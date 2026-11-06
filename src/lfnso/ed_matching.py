"""Conditional multi-zone ED position matching with occupancy-weighted Born factors.

No dynamical diffraction, detector intensity calibration, or gate pass is inferred.
Similarity registration uses the parent reciprocal grid, never a candidate's labels.
"""
import argparse, base64, csv, hashlib, html, json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment, minimize
from scipy.spatial import cKDTree, distance
from pymatgen.core import Lattice, Structure, Element
from pymatgen.io.cif import CifParser
from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS
from .model_generator import make_parent
from .pipeline import digest, write_json


def plane_reflections(parent_lattice, period, zone, gmax):
    """hkl refer to (a, period*b, c); zone is in the ORIGINAL parent direct basis."""
    z = np.asarray(zone, int)
    bounds = np.ceil(np.array(parent_lattice.abc)*np.array([1,period,1])*gmax/10).astype(int)+1
    # Solve the axis with a nonzero zone coefficient, avoiding a full 3D mesh.
    coefficient = np.array([period*z[0], z[1], period*z[2]])
    solved = int(np.flatnonzero(coefficient)[-1])
    free = [i for i in range(3) if i != solved]
    aa, bb = np.meshgrid(np.arange(-bounds[free[0]],bounds[free[0]]+1),np.arange(-bounds[free[1]],bounds[free[1]]+1),indexing='ij')
    hk = np.zeros((aa.size,3),int); hk[:,free[0]]=aa.ravel(); hk[:,free[1]]=bb.ravel()
    numerator = -(hk@coefficient)
    mask = numerator % coefficient[solved] == 0
    hk=hk[mask]; hk[:,solved]=numerator[mask]//coefficient[solved]
    parent_hkl=hk/np.array([1,period,1])
    g=parent_hkl@parent_lattice.reciprocal_lattice_crystallographic.matrix*10
    norms=np.linalg.norm(g,axis=1); mask=(norms>1e-8)&(norms<=gmax+1e-9)
    hk=hk[mask];g=g[mask];parent_hkl=parent_hkl[mask]
    normal=z@parent_lattice.matrix;normal=normal/np.linalg.norm(normal)
    seed=np.eye(3)[np.argmin(abs(normal))]
    e1=np.cross(normal,seed);e1/=np.linalg.norm(e1);e2=np.cross(normal,e1)
    xy=g@np.array([e1,e2]).T
    return hk,parent_hkl,g,xy


def electron_factor(symbol, g_nm_inv):
    """Neutral-atom Mott-Bethe factor in Å, stable at small g; no ionic/DW factors."""
    coeff=np.asarray(ATOMIC_SCATTERING_PARAMS[symbol])
    s2=(np.asarray(g_nm_inv)/20)**2
    return .023934*41.78214*np.sum(coeff[:,0]*np.exp(-s2[:,None]*coeff[:,1]),axis=1)


def electron_intensity(structure, hk, g_norm):
    fs={element.symbol:electron_factor(element.symbol,g_norm) for element in structure.elements}
    intensity=[]
    for indices in np.array_split(np.arange(len(hk)),max(1,int(np.ceil(len(hk)/512)))):
        phases=np.exp(2j*np.pi*(hk[indices]@structure.frac_coords.T))
        total=np.zeros(len(indices),complex)
        for i,site in enumerate(structure):
            # TEMCalculator.cell_scattering_factors ignores site occupancies; sum explicitly.
            for element,occupancy in site.species.items():total+=occupancy*fs[element.symbol][indices]*phases[:,i]
        intensity.extend(abs(total)**2)
    return np.asarray(intensity)


def transform(points, parameters, mirror=1):
    angle,scale,dy,dx=parameters
    rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    return (np.asarray(points)*np.array([1,mirror]))@rotation.T*scale+np.array([dy,dx])


def register_parent(observed, reference, scale_bounds, offset_bound, fit_count=32):
    """Robust trimmed positional objective; report alternatives, not unique indexing."""
    sample=np.asarray(observed[:fit_count]);keep=max(4,int(np.ceil(.65*len(sample))))
    tree=cKDTree(reference)
    def loss(p,mirror):
        angle,scale,dy,dx=p
        r=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        inverse=((sample-[dy,dx])/scale)@r*np.array([1,mirror])
        d=tree.query(inverse)[0]*scale
        return float(np.mean(np.partition(d*d,keep-1)[:keep]))
    solutions=[]
    for mirror in [1,-1]:
        coarse=[]
        for theta in np.deg2rad(np.arange(0,180,2)):
            for scale in np.linspace(*scale_bounds,11):
                p=[theta,scale,0,0];coarse.append((loss(p,mirror),p))
        starts=[]
        for _,p in sorted(coarse,key=lambda v:v[0]):
            if all(abs(np.angle(np.exp(2j*(p[0]-s[0]))))/2>np.deg2rad(8) for s in starts):starts.append(p)
            if len(starts)==4:break
        for p in starts:
            bounds=[(p[0]-np.deg2rad(5),p[0]+np.deg2rad(5)),tuple(scale_bounds),(-offset_bound,offset_bound),(-offset_bound,offset_bound)]
            fit=minimize(loss,p,args=(mirror,),method='Powell',bounds=bounds,options={'maxiter':150,'xtol':1e-7,'ftol':1e-8})
            solutions.append({'parameters':fit.x.tolist(),'mirror':mirror,'trimmed_RMS_nm_inv':float(np.sqrt(fit.fun)),'optimizer_success':bool(fit.success)})
    solutions.sort(key=lambda r:r['trimmed_RMS_nm_inv'])
    return {**solutions[0],'alternative_solutions':solutions[1:], 'fit_observed_count':len(sample),'trimmed_count':keep,
            'scale_definition':'g_in_metadata_units = scale * rotated_g_from_crystal + offset; not an applied metadata correction',
            'scale_bounds':scale_bounds,'center_offset_bound_nm_inv':offset_bound,'fit_source':'strongest/earliest detected positions on parent C-centered grid; labels/intensities of candidates unused'}


def match_one_to_one(observed,predicted,tolerance):
    """Maximum cardinality within tolerance, then minimum distance; every spot used once."""
    observed=np.asarray(observed);predicted=np.asarray(predicted)
    n,m=len(observed),len(predicted)
    if not n or not m:return {'matches':0,'observed_count':n,'predicted_count':m,'recall':0.,'precision':0.,'f1':0.,'matched_RMS_nm_inv':None,'pairs':[]}
    d=distance.cdist(observed,predicted)
    # A missed match costs more than all valid residual costs together.
    penalty=n+1.;cost=np.full((n,m+n),penalty)
    cost[:,:m]=np.where(d<=tolerance,d/tolerance,penalty*2)
    ii,jj=linear_sum_assignment(cost);valid=(jj<m)&(cost[ii,jj]<penalty)
    ii=ii[valid];jj=jj[valid];k=len(ii);precision=k/m;recall=k/n
    return {'matches':k,'observed_count':n,'predicted_count':m,'recall':recall,'precision':precision,
            'f1':2*k/(n+m),'matched_RMS_nm_inv':float(np.sqrt(np.mean(d[ii,jj]**2))) if k else None,
            'pairs':[{'observed_index':int(i),'predicted_index':int(j),'distance_nm_inv':float(d[i,j])} for i,j in zip(ii,jj)]}


def in_view(points, dataset):
    points=np.asarray(points);r=np.linalg.norm(points,axis=1)
    px=points/dataset['pixel_scale']+dataset['center']
    y0,x0,y1,x1=dataset['roi']
    return (r>=dataset['g_range'][0])&(r<=dataset['g_range'][1])&(px[:,0]>=y0)&(px[:,0]<y1)&(px[:,1]>=x0)&(px[:,1]<x1)


def load_dataset(entry,base):
    p=(base/entry['peaks']).resolve();record=json.loads(p.read_text())
    points=np.array([v['g_yx_nm_inv'] for v in record['peaks']]);scale=entry['pixel_scale_nm_inv']
    dataset={'id':entry['id'],'pixel_scale':scale,'center':np.array(record['center_yx_px']),'roi':entry['roi_yxyx'],
             'g_range':entry.get('g_range_nm_inv',[.6,12]),'raw':entry['raw'],'peak_source':str(p),'peak_source_sha256':digest(p),
             'raw_sha256':digest(Path(entry['raw'])),'nominal_zone':entry['nominal_zone'],'kind':entry['kind'],'entry':entry}
    retained=[];indices=[]
    for i,point in enumerate(points):
        if not in_view([point],dataset)[0]:continue
        if retained and min(np.linalg.norm(np.array(retained)-point,axis=1))<2*scale:continue
        retained.append(point);indices.append(i)
    dataset['observed']=np.array(retained);dataset['source_peak_indices']=indices
    dataset['position_tolerances_nm_inv']=entry['position_tolerances_nm_inv']
    dataset['filter_note']='g/ROI selection, then greedy centroid duplicate removal within 2 pixels in original detection order; not independent observations'
    return dataset


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from PIL import Image
    from .data_io import read_dm
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    write_json(out/'config.json',cfg)
    datasets=[load_dataset(e,cp.parent) for e in cfg['datasets']]
    source=(cp.parent/cfg['candidate_summary']).resolve();search=json.loads(source.read_text())
    template=Structure.from_file(cfg['template']);models=[]
    for row in search['rows']:
        folder=source.parent/row['id'];s=CifParser(folder/'model.cif',frac_tolerance=0).parse_structures(primitive=False)[0]
        models.append({'id':row['id'],'metric':row['metric_name'],'period':row['period_b'],'structure':s,'cif_sha256':digest(folder/'model.cif'),
                       'xrd_R2_full_superstructure':row['scores']['train']['R2_vs_linear_background'],'role':'trial'})
    for name,metric in search['config']['metrics'].items():
        models.append({'id':name+'_ordered_control','metric':name,'period':1,'structure':make_parent(template,[8,9],metric),'role':'control'})
    registrations={}; rows=[]; patterns={}; reference_cache={}
    for dataset in datasets:
        for metric_name,metric in search['config']['metrics'].items():
            lattice=Lattice.from_parameters(*metric)
            for zone in dataset['entry']['zone_hypotheses']:
                zone_key=','.join(map(str,zone));key=dataset['id']+'__'+metric_name+'__'+zone_key
                _,hp,_,xy=plane_reflections(lattice,1,zone,14)
                reference=xy[np.sum(hp[:,:2],axis=1).astype(int)%2==0]
                registration=register_parent(dataset['observed'],reference,cfg['registration_scale_bounds'],3*dataset['pixel_scale'])
                registrations[key]=registration;reference_cache[key]=reference
                print('registration',key,'scale',round(registration['parameters'][1],5),'rms',round(registration['trimmed_RMS_nm_inv'],5),flush=True)
                for model in models:
                    if model['metric']!=metric_name:continue
                    hk,hp,g,xy=plane_reflections(lattice,model['period'],zone,14)
                    intensity=electron_intensity(model['structure'],hk,np.linalg.norm(g,axis=1))
                    pred=transform(xy,registration['parameters'],registration['mirror']);visible=in_view(pred,dataset)
                    norm=max(float(intensity[visible].max()) if visible.any() else 0,1e-30)
                    relative=intensity/norm*100
                    pattern_key=key+'__'+model['id']
                    patterns[pattern_key]={'hkl_supercell':hk.tolist(),'hkl_parent':hp.tolist(),'g_projected_yx_nm_inv':pred.tolist(),'relative_kinematic_intensity':relative.tolist(),'inside_comparison_roi':visible.tolist()}
                    scores=[]
                    for threshold in cfg['intensity_thresholds_percent']:
                        selected=visible&(relative>=threshold);actual=pred[selected]
                        for tolerance in dataset['position_tolerances_nm_inv']:
                            score=match_one_to_one(dataset['observed'],actual,tolerance)
                            prediction_ids=np.flatnonzero(selected)
                            for pair in score['pairs']:
                                pair['source_peak_index']=dataset['source_peak_indices'][pair['observed_index']]
                                pair['pattern_reflection_index']=int(prediction_ids[pair['predicted_index']])
                            scores.append({'threshold_percent':threshold,'tolerance_nm_inv':tolerance,**score})
                    standard=next(s for s in scores if s['threshold_percent']==1 and s['tolerance_nm_inv']==dataset['position_tolerances_nm_inv'][1])
                    nulls=[]
                    for extra in cfg['orientation_null_degrees']:
                        p=registration['parameters'].copy();p[0]+=np.deg2rad(extra)
                        pp=transform(xy,p,registration['mirror']);vis=in_view(pp,dataset)
                        maximum=max(float(intensity[vis].max()) if vis.any() else 0,1e-30)
                        selected=vis&(intensity>=maximum*.01)
                        ns=match_one_to_one(dataset['observed'],pp[selected],standard['tolerance_nm_inv'])
                        nulls.append({'extra_rotation_degrees':extra,'f1':ns['f1'],'recall':ns['recall'],'precision':ns['precision']})
                    row={'dataset':dataset['id'],'model':model['id'],'role':model['role'],'metric':metric_name,'period_b':model['period'],
                         'zone_hypothesis':zone,'nominal_zone':dataset['nominal_zone'],'registration_key':key,'pattern_key':pattern_key,
                         'scores':scores,'standard_score':{k:v for k,v in standard.items() if k!='pairs'},'rotated_controls':nulls,
                         'xrd_R2_full_superstructure':model.get('xrd_R2_full_superstructure'),
                         'gate2':'not established; calibration/indexing/dynamical effects unresolved'}
                    rows.append(row)
    (out/'patterns').mkdir()
    for key,pattern in patterns.items():write_json(out/'patterns'/f'{key}.json',pattern)
    dataset_records=[{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in d.items() if k not in ['entry']} for d in datasets]
    summary={'config':cfg,'code_sha256':digest(Path(__file__)),'candidate_summary_sha256':digest(source),'template_sha256':digest(Path(cfg['template'])),
             'models':[{k:v for k,v in m.items() if k!='structure'} for m in models], 'datasets':dataset_records,'registrations':registrations,'rows':rows,
             'method':'Strict zero-order Laue-zone positions, neutral-atom Mott-Bethe factors with .67/.33 occupancies, one-to-one matches; all candidates share parent-based pose per dataset/metric/zone.',
             'limitations':['F1 is a descriptive position score, not a probability or crystallographic residual.',
                            'Detection thresholds/caps and merging affect observed/absent spots; predicted unmatched spots are not proof of systematic absences.',
                            'Kinematic intensities only choose threshold subsets. No fit to measured ED intensity or dynamical scattering performed.',
                            'Scale variation is a diagnostic fit, not a measured calibration correction. Mirror and zone alternatives are explicitly retained.',
                            'A larger reciprocal grid can increase recall by chance; precision and rotated same-model controls are retained. Controls are not statistical p-values.',
                            'TIFFs are display exports. Peaks, symmetry mates, images and fitted registration are not statistically independent.',
                            'Strict [110] ZOLZ of b-only supercells permits k_super=-period*h: this family cannot add fractional b satellites within that exact plane.',
                            'XRD gates remain unpassed. These exploratory ED diagnostics do not establish a final structure or authorize an unqualified STEM validation claim.'],
             'validated_shortlist':[],'stem_multislice':'not run'}
    write_json(out/'summary.json',summary)
    with (out/'scores.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['dataset','model','zone','threshold_percent','tolerance_nm_inv','matched','observed','predicted','recall','precision','F1','RMS_nm_inv'])
        for row in rows:
            for score in row['scores']:writer.writerow([row['dataset'],row['model'],str(row['zone_hypothesis']),score['threshold_percent'],score['tolerance_nm_inv'],score['matches'],score['observed_count'],score['predicted_count'],score['recall'],score['precision'],score['f1'],score['matched_RMS_nm_inv']])
    # Overlay the same XRD-best trial in every image, avoiding ED-cherry-picked model changes.
    selected_model=cfg['overlay_model'];overlays=[]
    for dataset in datasets:
        if dataset['kind']=='dm3':a,_=read_dm(Path(dataset['raw']));shown=np.log1p(np.maximum(a,0));lo,hi=np.percentile(shown,[1,99.8])
        else:shown=np.array(Image.open(dataset['raw']).convert('L'));lo,hi=0,255
        for zone in dataset['entry']['zone_hypotheses']:
            row=next(r for r in rows if r['model']==selected_model and r['dataset']==dataset['id'] and r['zone_hypothesis']==zone)
            pattern=patterns[row['pattern_key']];pred=np.array(pattern['g_projected_yx_nm_inv']);vis=np.array(pattern['inside_comparison_roi'])&(np.array(pattern['relative_kinematic_intensity'])>=1)
            pixels=pred[vis]/dataset['pixel_scale']+dataset['center'];obs=dataset['observed']/dataset['pixel_scale']+dataset['center']
            fig,ax=plt.subplots(figsize=(8,8));ax.imshow(shown,cmap='gray',vmin=lo,vmax=hi)
            ax.scatter(obs[:,1],obs[:,0],s=45,edgecolors='tomato',facecolors='none',lw=.8,label='Detected position')
            ax.scatter(pixels[:,1],pixels[:,0],s=20,marker='+',c='cyan',lw=.6,label='Predicted >=1% (kinematic)')
            y0,x0,y1,x1=dataset['roi'];ax.set(xlim=(x0,x1),ylim=(y1,y0),title=f"{dataset['id']} / zone hypothesis {zone}\n{selected_model}; F1={row['standard_score']['f1']:.3f}")
            ax.legend(fontsize=8);fig.tight_layout();name=dataset['id']+'__'+','.join(map(str,zone))+'.png';fig.savefig(out/name,dpi=140);plt.close(fig);overlays.append(name)
    page='<!doctype html><meta charset="utf-8"><title>LFNSO21-B ED matching</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto}img{max-width:100%}pre{white-space:pre-wrap}td,th{padding:6px;border-bottom:1px solid #ccc}</style><h1>실제 ED–후보 위치 정합 · 조건부 진단</h1><p>DM3 2장과 표시용 TIFF 2장을 비교했습니다. [101] 파일은 [101]/[10-1] 가설을 모두 보존했습니다. 회절 위치는 C-centered parent grid로 정합하고, 후보별 점유 배열에 따른 전자 산란인자 계산으로 반사 threshold를 적용했습니다. 다중 산란·두께·tilt를 피팅하지 않았으며 최종 gate 통과를 뜻하지 않습니다. 같은 XRD 최상 후보를 모든 overlay에 표시합니다.</p>'
    page+='<h2>등록 격자와 scale 진단</h2><table><tr><th>영상 / metric / zone</th><th>예측→metadata scale</th><th>trimmed RMS / nm⁻¹</th></tr>'
    for key,r in registrations.items():page+=f'<tr><td>{key}</td><td>{r["parameters"][1]:.5f}</td><td>{r["trimmed_RMS_nm_inv"]:.5f}</td></tr>'
    page+='</table>'
    for name in overlays:page+='<h2>'+html.escape(name)+'</h2><img src="data:image/png;base64,'+base64.b64encode((out/name).read_bytes()).decode()+'">'
    page+='<h2>수치·출처·모든 조건</h2><pre>'+html.escape(json.dumps(summary,indent=2,ensure_ascii=False))+'</pre>'
    (out/'report.html').write_text(page)
    (out/'source_snapshot').mkdir();(out/'source_snapshot'/'ed_matching.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({'output':str(out),'model_count':len(models),'comparison_count':len(rows)}),flush=True)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
