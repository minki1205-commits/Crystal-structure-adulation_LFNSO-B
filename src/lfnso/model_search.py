"""Small composition-preserving hypothesis library and first position screening."""
import argparse,csv,html,json,time,base64
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from pymatgen.core import Structure,Element
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .model_generator import make_parent,sb_orbits,generate,validate
from .pipeline import write_json,digest


def specifications():
    rows=[]
    for pattern in range(3):rows.append({'name':f'parent_pattern{pattern}','pattern':pattern,'ny':1,'nz':1,'family':'H0_reference'})
    for n,m in [(3,1),(4,1),(6,1),(7,2),(11,3),(12,3)]:
        for family,extra in [('exact',{}),('occupancy',{'exchange':.25}),('displacement',{'displacement_A':.02})]:
            rows.append({'name':f'{family}_{n}b_m{m}','pattern':0,'ny':n,'nz':1,'harmonic':m,'family':family,**extra})
    rows.append({'name':'stacking_AABB','pattern':0,'ny':1,'nz':4,'stacking':True,'family':'H2'})
    rows.append({'name':'coupled_4b_AABB','pattern':0,'ny':4,'nz':4,'harmonic':1,'exchange':.25,'stacking':True,'family':'H6'})
    rows.append({'name':'coupled_7b_AABB','pattern':0,'ny':7,'nz':4,'harmonic':2,'exchange':.25,'stacking':True,'family':'H6'})
    rows.append({'name':'integer_FeNi_3b_random','pattern':0,'ny':3,'nz':1,'integer_approx':True,'seed':0,'family':'integer_2over3_approximation'})
    for pattern in range(3):rows.append({'name':f'Sb_Fe_host_7b_pattern{pattern}','pattern':pattern,'ny':7,'nz':1,'harmonic':2,'host_exchange':1.,'family':'limited_Sb_Fe_host_control'})
    for name,extra in [('Li_TM_75',{'exchange':.75}),('Li_TM_100',{'exchange':1.}),('D008',{'displacement_A':.08}),('D015',{'displacement_A':.15})]:
        rows.append({'name':name+'_7b','pattern':0,'ny':7,'nz':1,'harmonic':2,'family':'amplitude_extension',**extra})
    return rows


def position_match(two_theta,intensities,observed,threshold):
    mask=np.asarray(intensities)>=threshold;pred=np.asarray(two_theta)[mask]
    if len(pred)<len(observed):return {'status':'too few predicted peaks above analysis threshold','matched':len(pred),'rms_deg':None}
    cost=abs(np.asarray(observed)[:,None]-pred[None,:]);ii,jj=linear_sum_assignment(cost)
    return {'status':'conditional geometric matching','matched':len(ii),'rms_deg':float(np.sqrt(np.mean(cost[ii,jj]**2))),
            'observed_deg':np.asarray(observed)[ii].tolist(),'predicted_deg':pred[jj].tolist(),'residual_deg':(pred[jj]-np.asarray(observed)[ii]).tolist(),
            'predicted_peak_count_above_threshold':len(pred),'threshold_relative_to_simulated_max_percent':threshold}


def proxy_ed(s,ny,nz,out):
    weights=np.array([sum(Element(e.symbol).Z*v for e,v in t.species.items()) for t in s])
    # Explicit approximation: constant Z weights, not electron scattering factors or dynamical ED.
    for zone in [[0,0,1],[1,0,1],[1,1,0]]:
        hk=[]
        for h in range(-2,3):
            for k in range(-3*ny,3*ny+1):
                for l in range(-2*nz,2*nz+1):
                    parent=np.array([h,k/ny,l/nz])
                    if np.any(parent) and abs(parent@zone)<1e-8:hk.append([h,k,l])
        hk=np.array(hk);g=hk@s.lattice.reciprocal_lattice_crystallographic.matrix*10
        mask=np.linalg.norm(g,axis=1)<=5;hk=hk[mask];g=g[mask]
        values=[]
        for chunk in np.array_split(hk,max(1,int(np.ceil(len(hk)/128)))):
            values.extend(abs(np.exp(2j*np.pi*(chunk@s.frac_coords.T))@weights)**2)
        vals=np.asarray(values);vals=vals/max(float(vals.max()),1)*100
        with (out/('zone_'+''.join(map(str,zone))+'.csv')).open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['h_super','k_super','l_super','h_parent','k_parent','l_parent','g_x_nm-1','g_y_nm-1','g_z_nm-1','Z_proxy_relative_intensity_NOT_electron_intensity'])
            for h,v,i in zip(hk,g,vals):w.writerow([*h,*(h/np.array([1,ny,nz])),*v,i])


def run(config_path):
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    write_json(out/'config.json',cfg);template=Structure.from_file(cfg['template']);reps=sb_orbits(template)
    parents=[make_parent(template,r,cfg['metric_parameters']) for r in reps]
    raw=json.loads((cp.parent/cfg['xrd_result']).resolve().read_text());fs=raw['fits']['LFNSO21-B']
    labels=['anchor_001_assumed','anchor_020_assumed','feature_1','feature_2','feature_3']
    observed=np.array([fs[k]['center_two_theta_deg'] for k in labels])
    rows=[];t0=time.monotonic();reference_pattern=None
    for spec in specifications():
        if cfg.get("candidate_ids") and spec["name"] not in cfg["candidate_ids"]:continue
        name=spec['name'];folder=out/name;folder.mkdir()
        args={k:v for k,v in spec.items() if k not in ['name','pattern','family']}
        s=generate(parents[spec['pattern']],**args);n=spec['ny']*spec['nz'];v=validate(s,2*n)
        for site in s:site.label=site.species_string
        s.to(filename=str(folder/'model.cif'),fmt='cif',symprec=None)
        record={'id':name,'specification':spec,'template_sha256':digest(Path(cfg['template'])),'Sb_parent_site_indices':list(reps[spec['pattern']]),
                'generation_note':'Li and O coordinates inherited from supplied Li2MnO3; Mn sublattice divided into Sb and Fe0.67Ni0.33 host sites. P1 is computational, not an established space group.',
                'transformation':np.diag([1,spec['ny'],spec['nz']]).tolist(),'site_count':len(s),'validation':v,
                'relaxation':'not performed','gate1':'not_evaluated','gate2':'not_evaluated','stem':'not_run'}
        pattern=XRDCalculator(wavelength=cfg['wavelength_A']).get_pattern(s,two_theta_range=(17,36))
        write_json(folder/'xrd.json',{'two_theta_deg':pattern.x.tolist(),'relative_intensity':pattern.y.tolist(),'wavelength_A':cfg['wavelength_A'],'wavelength_status':'provisional Cu Kalpha','interpretation':'ideal periodic average structure; no instrument profile or disorder diffuse scattering'})
        record['position_screen']={str(c):position_match(pattern.x,pattern.y,observed,c) for c in [.1,.5,1.]}
        if name=='parent_pattern0':reference_pattern=(pattern.x.copy(),pattern.y.copy())
        if spec['family']=='exact':
            xx,yy=reference_pattern
            record['repeat_control_pass']=bool(len(xx)==len(pattern.x) and np.allclose(xx,pattern.x,atol=1e-6) and np.allclose(yy,pattern.y,atol=1e-5))
            if not record['repeat_control_pass']:raise AssertionError('Exact repetition changed powder diffraction')
        proxy_ed(s,spec['ny'],spec['nz'],folder)
        write_json(folder/'metadata.json',record);rows.append(record)
        print(name,'sites',len(s),'geometry',v['gross_geometry_pass'],'CN6',round(v['sixfold_O_coordination_fraction'],3),flush=True)
    summary={'candidate_count':len(rows),'elapsed_seconds':time.monotonic()-t0,'Sb_symmetry_distinct_patterns':[list(r) for r in reps],
             'experimental_XRD_source':str((cp.parent/cfg['xrd_result']).resolve()),'observed_features':labels,
             'candidates':[{k:v for k,v in r.items() if k!='validation'}|{'validation_summary':{k:v for k,v in r['validation'].items() if k!='coordination'}} for r in rows],
             'limitations':['Partial occupancy diffraction is an average model, not a specific integer atomic configuration.',
                            'All metric parameters fixed to the master-prompt starting cell; position mismatches cannot distinguish site ordering by themselves.',
                            'Only Sb pattern0 modulation family explored; other base patterns and symmetry-reduced supercell enumeration remain.',
                            'Proxy ED Z weights do not predict experimental electron diffraction intensities.',
                            'XRD thresholds are analysis sensitivity settings, not experimental detection limits. Dense supercells can match more peaks by chance.',
                            'Oxygen coordinates unrelaxed; ionic charge and local geometry remain conditional.',
                            'No final model ranking, gate pass or STEM simulation is claimed.']}
    write_json(out/'summary.json',summary)
    table=[]
    for r in rows:
        score=r['position_screen']['0.5'];table.append([r['id'],r['site_count'],r['validation']['exact_decimal_composition'],r['validation']['sixfold_O_coordination_fraction'],score['rms_deg'],r.get('repeat_control_pass','')])
    with (out/'screening.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['id','sites','exact_Li4Fe0.67Ni0.33SbO6','CN6_fraction_at_2p7A','position_RMS_deg_at_0p5percent_threshold','repeat_control_pass']);w.writerows(table)
    fig,ax=plt.subplots(figsize=(11,5));x=np.arange(len(rows));ys=[r['position_screen']['0.5']['rms_deg'] for r in rows]
    ax.bar(x,[v if v is not None else np.nan for v in ys],color='steelblue');ax.set_xticks(x,[r['id'] for r in rows],rotation=80,fontsize=7)
    ax.set(ylabel='Conditional position RMS (degrees)',title='First hypothesis library — not an atomic structure ranking')
    fig.tight_layout();fig.savefig(out/'position_screen.png',dpi=150);plt.close(fig)
    page='<!doctype html><meta charset="utf-8"><title>LFNSO hypothesis library</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto}pre{white-space:pre-wrap}img{max-width:100%}</style><h1>LFNSO 가설 구조 생성 및 첫 위치 screening</h1><p>사용자 조성과 Fe-site Ni 치환 가정을 반영했습니다. 원자 배열·산화수·적층 구조는 검증 전 가설입니다. 평균 부분 점유 모델과 정수 배열 근사를 구분합니다.</p>'
    page+='<img src="data:image/png;base64,'+base64.b64encode((out/'position_screen.png').read_bytes()).decode()+'"><pre>'+html.escape(json.dumps(summary,indent=2,ensure_ascii=False))+'</pre>'
    (out/'report.html').write_text(page);return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
