"""Bounded affine diagnostic of ED metric/registration confounding, not calibration."""
import argparse,base64,csv,html,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from pymatgen.core import Lattice
from .ed_matching import plane_reflections,transform,match_one_to_one,in_view
from .pipeline import write_json,digest


def apply_correction(points,correction):
    return (np.asarray(points)-correction['base_offset'])@np.array(correction['matrix']).T+correction['total_offset']


def fit_correction(observed,reference,base_offset,offset_bound,tolerance,robust_scale,component_bound=.04):
    sample=np.asarray(observed[:32]);base_offset=np.asarray(base_offset)
    p=np.r_[np.zeros(4),base_offset]
    lower=np.r_[np.full(4,-component_bound),np.full(2,-offset_bound)]
    upper=-lower
    def corrected(p,points):return (points-base_offset)@(np.eye(2)+p[:4].reshape(2,2)).T+p[4:]
    for _ in range(4):
        matches=match_one_to_one(sample,corrected(p,reference),tolerance)
        pairs=matches['pairs']
        if len(pairs)<6:raise ValueError('Too few parent correspondences for affine diagnostic')
        ii=np.array([v['observed_index'] for v in pairs]);jj=np.array([v['predicted_index'] for v in pairs])
        if np.linalg.matrix_rank(np.c_[reference[jj],np.ones(len(jj))])<3:raise ValueError('Parent correspondences are collinear')
        fit=least_squares(lambda pp:(corrected(pp,reference[jj])-sample[ii]).ravel(),p,bounds=(lower,upper),loss='soft_l1',f_scale=robust_scale,max_nfev=300)
        p=fit.x
    matrix=np.eye(2)+p[:4].reshape(2,2)
    final=match_one_to_one(sample,corrected(p,reference),tolerance)
    residual=np.array([v['distance_nm_inv'] for v in final['pairs']]);count=max(4,int(np.ceil(.65*len(sample))))
    # Fixed-count trimmed nearest residual is directly comparable to initial registration loss.
    from scipy.spatial import cKDTree
    d=cKDTree(corrected(p,reference)).query(sample)[0]
    return {'matrix':matrix.tolist(),'total_offset':p[4:].tolist(),'base_offset':base_offset.tolist(),
            'incremental_singular_values':np.linalg.svd(matrix,compute_uv=False).tolist(),
            'component_bound':component_bound,'component_bound_hit':bool(np.any(abs(p[:4])>component_bound-1e-5)),
            'total_offset_bound_nm_inv':offset_bound,'offset_bound_hit':bool(np.any(abs(p[4:])>offset_bound-1e-5)),
            'parent_correspondence_count':len(final['pairs']),'parent_correspondence_RMS_nm_inv':float(np.sqrt(np.mean(residual**2))),
            'trimmed_RMS_nm_inv':float(np.sqrt(np.mean(np.partition(d*d,count-1)[:count]))),
            'fit_source':'same first 32 detections, robust assignment to C-centered parent grid; candidate occupations unused',
            'interpretation':'Sensitivity only: anisotropic calibration, different local metric, tilt and mistaken correspondence are confounded; not a calibrated diffraction fit.'}


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from PIL import Image
    from .data_io import read_dm
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());source=(cp.parent/cfg['source_summary']).resolve();s=json.loads(source.read_text())
    out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    write_json(out/'config.json',cfg)
    model_metrics={m['metric']:None for m in s['models']}
    search_path=(cp.parent/s['config']['candidate_summary']).resolve();search=json.loads(search_path.read_text())
    model_metrics=search['config']['metrics'];corrections={};rows=[];data_by_id={d['id']:d for d in s['datasets']}
    for key,reg in s['registrations'].items():
        dataset_id,metric,zone_text=key.split('__');dataset=data_by_id[dataset_id];zone=list(map(int,zone_text.split(',')))
        _,hp,_,xy=plane_reflections(Lattice.from_parameters(*model_metrics[metric]),1,zone,14)
        reference=transform(xy[np.sum(hp[:,:2],axis=1).astype(int)%2==0],reg['parameters'],reg['mirror'])
        correction=fit_correction(dataset['observed'],reference,reg['parameters'][2:],3*dataset['pixel_scale'],.35 if dataset['kind']=='dm3' else .48,.05 if dataset['kind']=='dm3' else .16,cfg['affine_component_bound'])
        corrections[key]=correction
        print(key,'RMS',round(reg['trimmed_RMS_nm_inv'],4),'->',round(correction['trimmed_RMS_nm_inv'],4),'stretch',correction['incremental_singular_values'],flush=True)
    for row in s['rows']:
        dataset=data_by_id[row['dataset']];pattern=json.loads((source.parent/'patterns'/(row['pattern_key']+'.json')).read_text())
        corrected=apply_correction(pattern['g_projected_yx_nm_inv'],corrections[row['registration_key']]);visible=in_view(corrected,dataset)
        intensity=np.array(pattern['relative_kinematic_intensity']);norm=max(float(intensity[visible].max()) if visible.any() else 0,1e-30);relative=100*intensity/norm
        scores=[]
        for threshold in s['config']['intensity_thresholds_percent']:
            selected=visible&(relative>=threshold);ids=np.flatnonzero(selected)
            for tolerance in dataset['position_tolerances_nm_inv']:
                match=match_one_to_one(dataset['observed'],corrected[selected],tolerance)
                for pair in match['pairs']:
                    pair['source_peak_index']=dataset['source_peak_indices'][pair['observed_index']]
                    pair['pattern_reflection_index']=int(ids[pair['predicted_index']])
                scores.append({'threshold_percent':threshold,'tolerance_nm_inv':tolerance,**match})
        standard=next(v for v in scores if v['threshold_percent']==1 and v['tolerance_nm_inv']==dataset['position_tolerances_nm_inv'][1])
        rows.append({k:v for k,v in row.items() if k not in ['scores','standard_score','rotated_controls']}|{'similarity_standard_score':row['standard_score'],'affine_diagnostic_scores':scores,'affine_standard_score':{k:v for k,v in standard.items() if k!='pairs'}})
    summary={'source_summary_sha256':digest(source),'code_sha256':digest(Path(__file__)),'ed_matching_code_sha256':digest(Path(__file__).with_name('ed_matching.py')),
             'config':cfg,'corrections':corrections,'rows':rows,'validated_shortlist':[],
             'limitations':['This is an adaptive diagnostic after inspecting similarity-registration residuals, not an independent validation.',
                            'An added affine map can absorb real metric differences or detector distortion. It cannot decide their cause.',
                            '±4% bounds apply to each incremental matrix component, not a statistical calibration error; singular stretches are reported separately.',
                            'All candidates share the same parent-derived correction for each image/metric/zone. CIF coordinates and source images are unchanged.',
                            'Kinematic intensity threshold is arbitrary; no experimental ED intensity fit. Earlier rotated controls pertain only to similarity mode.',
                            'No joint XRD/ED winner or STEM gate pass is established.']}
    write_json(out/'summary.json',summary)
    with (out/'scores.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['dataset','model','zone','threshold_percent','tolerance_nm_inv','matches','observed','predicted','recall','precision','F1','RMS_nm_inv'])
        for row in rows:
            for v in row['affine_diagnostic_scores']:writer.writerow([row['dataset'],row['model'],str(row['zone_hypothesis']),v['threshold_percent'],v['tolerance_nm_inv'],v['matches'],v['observed_count'],v['predicted_count'],v['recall'],v['precision'],v['f1'],v['matched_RMS_nm_inv']])
    page='<!doctype html><meta charset="utf-8"><title>ED metric sensitivity</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto}img{max-width:100%}td,th{padding:6px;border-bottom:1px solid #ddd}pre{white-space:pre-wrap}</style><h1>ED 위치 정합 — 격자·영상 왜곡 민감도</h1><p>회전+등방 배율 정합 후 방향별 변형을 제한적으로 허용한 진단입니다. 실제 calibration이나 정련된 격자가 아닙니다. 모든 후보에 동일한 parent-derived 보정을 적용하고 원본·CIF는 변경하지 않았습니다.</p>'
    page+='<table><tr><th>영상 / metric / zone</th><th>기본 RMS</th><th>진단 RMS</th><th>추가 singular stretches</th><th>행렬 경계 도달</th></tr>'
    for key,c in corrections.items():page+=f'<tr><td>{key}</td><td>{s["registrations"][key]["trimmed_RMS_nm_inv"]:.4f}</td><td>{c["trimmed_RMS_nm_inv"]:.4f}</td><td>{c["incremental_singular_values"]}</td><td>{c["component_bound_hit"]}</td></tr>'
    page+='</table>'
    for dataset in s['datasets']:
        for row in rows:
            if row['dataset']!=dataset['id'] or row['model']!=s['config']['overlay_model']:continue
            if dataset['kind']=='dm3':a,_=read_dm(Path(dataset['raw']));a=np.log1p(np.maximum(a,0));lo,hi=np.percentile(a,[1,99.8])
            else:a=np.array(Image.open(dataset['raw']).convert('L'));lo,hi=0,255
            pattern=json.loads((source.parent/'patterns'/(row['pattern_key']+'.json')).read_text());pred=apply_correction(pattern['g_projected_yx_nm_inv'],corrections[row['registration_key']]);vis=in_view(pred,dataset);ii=np.array(pattern['relative_kinematic_intensity']);vis &= ii>=max(float(ii[vis].max()),1e-30)*.01
            px=pred[vis]/dataset['pixel_scale']+dataset['center'];obs=np.array(dataset['observed'])/dataset['pixel_scale']+dataset['center']
            fig,ax=plt.subplots(figsize=(8,8));ax.imshow(a,cmap='gray',vmin=lo,vmax=hi);ax.scatter(obs[:,1],obs[:,0],s=45,edgecolors='tomato',facecolors='none',lw=.8,label='Detected position');ax.scatter(px[:,1],px[:,0],marker='+',s=20,c='cyan',lw=.7,label='Predicted >=1%; affine diagnostic')
            y0,x0,y1,x1=dataset['roi'];ax.set(xlim=(x0,x1),ylim=(y1,y0),title=f"{dataset['id']} / {row['zone_hypothesis']} / affine diagnostic\n{row['model']}; F1={row['affine_standard_score']['f1']:.3f}");ax.legend(fontsize=8);fig.tight_layout()
            name=dataset['id']+'__'+','.join(map(str,row['zone_hypothesis']))+'.png';fig.savefig(out/name,dpi=120);plt.close(fig)
            page+='<h2>'+name+'</h2><img src="data:image/png;base64,'+base64.b64encode((out/name).read_bytes()).decode()+'">'
    page+='<pre>'+html.escape(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre>';(out/'report.html').write_text(page)
    (out/'ed_metric_diagnostic_snapshot.py').write_bytes(Path(__file__).read_bytes())
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
