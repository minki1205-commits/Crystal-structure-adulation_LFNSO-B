"""Descriptive multi-objective screening; no automatic scientific gate pass."""
import argparse,json,math,csv,html
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from .pipeline import write_json


def nondominated(rows):
    result=[]
    for i,row in enumerate(rows):
        a=np.array(row['objectives'])
        if not any(np.all(np.array(other['objectives'])<=a+np.array([1e-6,1e-9,0])) and np.any(np.array(other['objectives'])<a-np.array([1e-6,1e-9,0])) for j,other in enumerate(rows) if j!=i):result.append(row['id'])
    return result


def run(model_dir,experimental_csv):
    root=Path(model_dir);destination=root/'assessment.json'
    if destination.exists():raise FileExistsError(destination)
    summary=json.loads((root/'summary.json').read_text());a=np.loadtxt(experimental_csv,delimiter=',',skiprows=1)
    mask=(a[:,0]>=20.4)&(a[:,0]<=34.6);x,y=a[mask].T
    rows=[]
    for r in summary['candidates']:
        p=json.loads((root/r['id']/'xrd.json').read_text());t=np.array(p['two_theta_deg']);inten=np.array(p['relative_intensity'])
        scores=[]
        for width in [.15,.25,.4]:
            profile=(np.exp(-4*np.log(2)*((x[:,None]-t[None,:])/width)**2)@inten)/width
            design=np.c_[profile,np.ones(len(x)),x-x.mean()]
            fit=lsq_linear(design,y,bounds=([0,-np.inf,-np.inf],[np.inf,np.inf,np.inf]))
            residual=design@fit.x-y
            scores.append({'gaussian_fwhm_deg':width,'normalized_RMS_by_observed_range':float(np.sqrt(np.mean(residual**2))/np.ptp(y)),
                           'scale_background_slope':fit.x.tolist()})
        spec=r['specification'];v=r['validation_summary'];active=bool(spec.get('exchange') or spec.get('host_exchange') or spec.get('displacement_A') or spec.get('integer_approx'))
        repeat=spec['ny']//math.gcd(spec.get('harmonic',1),spec['ny']) if active else 1
        complexity=repeat*(4 if spec.get('stacking') else 1)
        error=r['position_screen']['0.5']['rms_deg']
        eligible=v['exact_decimal_composition'] and v['gross_geometry_pass'] and v['sixfold_O_coordination_fraction']==1 and error is not None and any(t['scale_background_slope'][0]>1e-8 for t in scores)
        row={'profile_signal_supported':any(t['scale_background_slope'][0]>1e-8 for t in scores),'id':r['id'],'profile_scores':scores,'physical_repeat_cells_surrogate':complexity,
             'eligible_for_descriptive_Pareto':eligible,'objectives':[error if error is not None else 999.,min(s['normalized_RMS_by_observed_range'] for s in scores),complexity]}
        rows.append(row)
    pareto=nondominated([r for r in rows if r['eligible_for_descriptive_Pareto']])
    result={'rows':rows,'nondominated_ids':pareto,'objectives':['five-feature position RMS','best profile RMS among three explicitly trialled widths','physical repetition cell count surrogate'],
            'profile_window_deg':[20.4,34.6],'experimental_csv':str(Path(experimental_csv).resolve()),
            'interpretation':'Exploratory Pareto shortlist, not a model probability or validated atomic structure ranking.',
            'limitations':['Uniform Gaussian width is a nuisance sensitivity trial, not instrument calibration.',
                           'Profile scale and linear background fitted; no peak-specific shifts or intensity adjustment.',
                           'Position-only scores favor denser patterns; profile residual and repetition cost are separate checks.',
                           'Partial occupancy average structures are not local configurations; charge scenarios remain unresolved.',
                           'Coordination filter at 2.7 Å is a conservative screening rule, not proof of chemical stability.',
                           'Equivalent exact-repeat controls remain tied; degeneracy is not independent model evidence.']}
    write_json(destination,result)
    with (root/'pareto_scores.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['id','position_RMS_deg','profile_normalized_RMS','repeat_cell_surrogate','geometry_composition_filter','Pareto_nondominated'])
        for r in rows:w.writerow([r['id'],*r['objectives'],r['eligible_for_descriptive_Pareto'],r['id'] in pareto])
    page=(root/'report.html').read_text()+'<h2>Exploratory multi-objective assessment</h2><pre>'+html.escape(json.dumps(result,indent=2,ensure_ascii=False))+'</pre>'
    (root/'report_assessed.html').write_text(page)
    print('Pareto shortlist:',pareto)
    for r in sorted(rows,key=lambda r:r['objectives'][1])[:8]:print(r['id'],r['objectives'],r['eligible_for_descriptive_Pareto'])
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('model_dir');p.add_argument('experimental_csv');a=p.parse_args();run(a.model_dir,a.experimental_csv)
