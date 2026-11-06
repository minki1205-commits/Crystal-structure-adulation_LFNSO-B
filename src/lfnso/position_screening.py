"""Position-only hypothesis checks with explicit provenance and conditional indexing.

No XRD intensity fits, chemical structure rankings or gate passes are inferred.
"""
import argparse,base64,csv,html,json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter,maximum_filter
from scipy.spatial import cKDTree
from scipy.optimize import minimize_scalar
from pymatgen.core import Lattice
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .pipeline import digest,write_json


def predicted_d(lattice, hkl, delta):
    hkl=np.asarray(hkl,float).copy();hkl[:,1]+=delta
    return np.array([lattice.d_hkl(h) for h in hkl])


def xrd_screen(lattice, hkl, observed, deltas, halfwidth):
    observed=np.asarray(observed,float)
    rows=[]
    for delta in deltas:
        pred=predicted_d(lattice,hkl,delta);err=pred-observed
        rows.append({'delta_bstar':float(delta),'predicted_d_A':pred.tolist(),'residual_d_A':err.tolist(),
                     'rms_d_A':float(np.sqrt(np.mean(err**2))),
                     'within_manual_reading_intervals':int(np.sum(np.abs(err)<=halfwidth))})
    opt=minimize_scalar(lambda d: np.mean((predicted_d(lattice,hkl,d)-observed)**2),bounds=(0,.5),method='bounded')
    return {'discrete_hypotheses':sorted(rows,key=lambda r:r['rms_d_A']),
            'free_delta':float(opt.x),'free_rms_d_A':float(np.sqrt(opt.fun)),
            'basis_assumptions':'fixed metric, supplied hkl families, first-order +b* satellites',
            'uncertainty_note':'Manual halfwidth is a visual reading bound, not Gaussian sigma. No confidence level assigned.',
            'interpretation':'Ranks conditional position hypotheses only; not atomic structure candidates.'}


def empirical_satellites(peak_data, seed_indices):
    g=np.array([p['g_yx_nm_inv'] for p in peak_data['peaks']])
    amp=np.array([p['background_subtracted_peak'] for p in peak_data['peaks']])
    p1,p2=g[seed_indices];basis=np.array([(p1-p2)/2,(p1+p2)/2])
    for _ in range(3):
        frac=g@np.linalg.inv(basis);h=np.round(frac[:,0]);k=2*np.round((frac[:,1]-h)/2)+h
        parent=(abs(frac[:,0]-h)<.08)&(abs(frac[:,1]-k)<.1)&(np.linalg.norm(g,axis=1)>.2)
        if np.linalg.matrix_rank(np.c_[h,k][parent])!=2:raise ValueError('Insufficient independent parent reflections')
        basis=np.linalg.lstsq(np.c_[h,k][parent],g[parent],rcond=None)[0]
    frac=g@np.linalg.inv(basis);h=np.round(frac[:,0]);k=2*np.round((frac[:,1]-h)/2)+h
    residual=frac-np.c_[h,k]
    # Broad candidate range includes 1/6, 1/4 and 1/3. Report selection explicitly.
    satellite=(abs(residual[:,0])<.05)&(abs(residual[:,1])>.15)&(abs(residual[:,1])<.45)&(amp>=50)
    if satellite.sum()<3:raise ValueError('Insufficient satellite-like peaks')
    signed=residual[satellite]*np.sign(residual[satellite,1])[:,None]
    delta=float(np.median(signed[:,1]));q=np.median(signed,axis=0)@basis
    records=[]
    for i in np.flatnonzero(satellite):
        records.append({'peak_index':int(i),'empirical_coordinates_uv':frac[i].tolist(),
                        'assigned_parent_uv':[float(h[i]),float(k[i])],
                        'residual_uv':residual[i].tolist(),'amplitude_detection_only':float(amp[i])})
    rng=np.random.default_rng(0)
    # Resampling estimates scatter only; pair symmetry and shared basis prevent independent statistical CI claims.
    medians=[float(np.median(rng.choice(signed[:,1],len(signed),replace=True))) for _ in range(1000)]
    grouped={}
    for lo,hi in [(0,4),(4,8),(8,100)]:
        mask=satellite&(np.abs(k)>=lo)&(np.abs(k)<hi)
        grouped[f'abs_parent_k_{lo}_{hi}']={'count':int(mask.sum()),'median_abs_delta':float(np.median(np.abs(residual[mask,1]))) if mask.any() else None}
    return {'basis_uv_rows_yx_nm_inv':basis.tolist(),'basis_norms_nm_inv':np.linalg.norm(basis,axis=1).tolist(),
            'basis_angle_deg':float(np.degrees(np.arccos(np.dot(*basis)/np.prod(np.linalg.norm(basis,axis=1))))),
            'parent_count':int(parent.sum()),'satellite_candidate_count':int(satellite.sum()),
            'median_delta_v':delta,'median_transverse_delta_u':float(np.median(signed[:,0])),
            'delta_v_robust_scatter':float(1.4826*np.median(abs(signed[:,1]-delta))),
            'bootstrap_median_2p5_97p5_descriptive_only':np.percentile(medians,[2.5,97.5]).tolist(),
            'q_image_yx_nm_inv':q.tolist(),'q_norm_nm_inv':float(np.linalg.norm(q)),
            'range_check':grouped,
            'selection':{'parent_u_tolerance':.08,'parent_v_tolerance':.1,'satellite_u_tolerance':.05,'satellite_v_range':[.15,.45],'min_detection_amplitude':50},
            'satellites':records,
            'limitations':['Basis was inferred from same image; indices are empirical, not validated crystallographic hkl.',
                          'Satellite selection is conditional on integer u and C-centered-like parent parity.',
                          'Selected opposite peaks are correlated, not independent measurements.',
                          'Scattering dynamics, overlapping domains and local distortions are not fitted.',
                          'Free delta is not evidence for a uniquely commensurate period.']}


def read_display_tiff(entry):
    path=Path(entry['path']);im=Image.open(path);a=np.asarray(im.convert('L'),float)
    y0,x0,y1,x1=entry['roi_yxyx'];xleft,xright=entry['scale_bar_x_edges'];row=entry['scale_bar_row']
    if not np.all(a[row,xleft:xright]>240):raise ValueError('Configured white scale bar not found')
    scale=entry['scale_bar_nm_inv']/(xright-xleft)
    r=gaussian_filter(a,.7)-gaussian_filter(a,5)
    mask=np.zeros_like(a,bool);mask[y0:y1,x0:x1]=True
    threshold=max(4,float(np.std(r[mask]))*1.2)
    pk=(r==maximum_filter(r,size=5))&(r>threshold)&mask
    y,x=np.where(pk);order=np.argsort(r[y,x])[::-1][:160]
    pts=np.c_[y[order],x[order]].astype(float)
    # Refine display centroids; intensity values are not raw electron counts.
    for i,(py,px) in enumerate(pts.copy()):
        py,px=int(py),int(px);patch=np.maximum(r[py-2:py+3,px-2:px+3],0)
        gy,gx=np.mgrid[py-2:py+3,px-2:px+3]
        if patch.sum():pts[i]=[(patch*gy).sum()/patch.sum(),(patch*gx).sum()/patch.sum()]
    mids=np.array([(pts[i]+pts[j])/2 for i in range(len(pts)) for j in range(i+1,len(pts))])
    seed=np.asarray(entry['center_seed_yx']);mids=mids[np.linalg.norm(mids-seed,axis=1)<8]
    if not len(mids):raise ValueError('No geometric pairs near the manually seeded beam center')
    cells,count=np.unique(np.round(mids).astype(int),axis=0,return_counts=True)
    center=cells[np.argmax(count)].astype(float);tree=cKDTree(pts)
    for _ in range(3):
        distances,indices=tree.query(2*center-pts)
        pairs={(min(i,int(j)),max(i,int(j))) for i,j in enumerate(indices) if i!=j and distances[i]<3}
        mm=np.array([(pts[i]+pts[j])/2 for i,j in pairs]);center=np.median(mm,axis=0)
    vecs=(pts-center)*scale
    report={'path':str(path),'sha256':digest(path),'zone_label':entry['zone_label'],'mode':im.mode,'shape_yx':list(a.shape),
            'calibration':'manual scale-bar annotation; equal x/y scale assumed for unwarped export',
            'scale_nm_inv_per_pixel':scale,'scale_bar_length_px':xright-xleft,'scale_bar_nm_inv':entry['scale_bar_nm_inv'],
            'bar_endpoint_width_uncertainty_px':2,'scale_relative_reading_bound':2/(xright-xleft),
            'center_yx_px':center.tolist(),'center_seed_yx':seed.tolist(),'opposite_pairs':len(pairs),
            'pair_midpoint_rms_yx_px':np.sqrt(np.mean((mm-center)**2,axis=0)).tolist(),
            'roi_yxyx':entry['roi_yxyx'],'threshold':threshold,'peak_count':len(pts),
            'peaks':[{'y_px':float(p[0]),'x_px':float(p[1]),'g_yx_nm_inv':v.tolist()} for p,v in zip(pts,vecs)],
            'quantitative_limitations':['8-bit palette display export; intensity/width not a raw diffraction measurement.',
                                       'Scale bars and label areas were excluded using recorded ROI.',
                                       'Small satellites can blend at reduced resolution; no hkl or q assigned from these TIFFs.']}
    return a,report


def run(config_path):
    path=Path(config_path).resolve();cfg=json.loads(path.read_text());out=(path.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    write_json(out/'config.json',cfg)
    lattice=Lattice.from_parameters(*cfg['lattice_parameters']);test=cfg['satellite_test']
    xrd=xrd_screen(lattice,test['parent_hkl'],test['B_d_A'],test['delta_candidates'],cfg['xrd_reading_halfwidth_A'])
    # An independent monoclinic metric identity conditional on the proposed hkl correspondence.
    d001,d020=cfg['xrd_B_d_A'][:2]
    predicted021=1/np.sqrt(d001**-2+d020**-2)
    w=cfg['xrd_reading_halfwidth_A']
    identity={'assumed_001_d_A':d001,'assumed_020_d_A':d020,'predicted_021_d_A':float(predicted021),
              'reading_interval_predicted_021_A':[float(1/np.sqrt((d001+sign*w)**-2+(d020+sign*w)**-2)) for sign in [-1,1]],
              'B_feature_tested_as_021_A':3.08,
              'conclusion':'If these hkl assignments are retained, the 3.08 Å feature cannot be the same monoclinic (021) reflection. This does not exclude another hkl or a metric-distorted alternative structure.'}
    peak_path=(path.parent/cfg['ed_peak_json']).resolve();peak_data=json.loads(peak_path.read_text());ed=empirical_satellites(peak_data,cfg['ed_seed_peak_indices'])
    ed['source_peak_file']=str(peak_path);ed['source_peak_sha256']=digest(peak_path)
    # This dimensionless ratio is independent of an overall detector calibration scale.
    observed_ratio=ed['basis_norms_nm_inv'][0]/ed['basis_norms_nm_inv'][1]
    ed['zone_convention_check']={'observed_u_over_v':observed_ratio,
         'expected_101_g10minus1_over_bstar':lattice.b/lattice.d_hkl([1,0,-1]),
         'expected_10minus1_g101_over_bstar':lattice.b/lattice.d_hkl([1,0,1]),
         'note':'Test only if v=b* and u is the shortest corresponding in-plane axial vector. Monoclinic axis sign/setting or indexing needs confirmation; do not silently relabel.'}
    sensitivity=[]
    for scale in [.99,1.,1.01]:
        lp=cfg['lattice_parameters'];scaled=Lattice.from_parameters(*(list(np.array(lp[:3])*scale)+lp[3:]))
        check=xrd_screen(scaled,test['parent_hkl'],test['B_d_A'],test['delta_candidates'],w)
        sensitivity.append({'length_scale':scale,'free_delta':check['free_delta'],'rms_d_A':check['free_rms_d_A']})
    pred=predicted_d(lattice,test['parent_hkl'],ed['median_delta_v']);obs=np.asarray(test['B_d_A'])
    length_scale=float(np.dot(pred,obs)/np.dot(pred,pred))
    predicted_bstar=10/(lattice.b*length_scale)
    joint={'delta_from_empirical_ED':ed['median_delta_v'],'fitted_uniform_length_scale':length_scale,
           'predicted_B_d_A':(pred*length_scale).tolist(),'rms_d_A':float(np.sqrt(np.mean((pred*length_scale-obs)**2))),
           'held_out_assumed_001_020_d_A':[float(lattice.d_hkl(h)*length_scale) for h in [[0,0,1],[0,2,0]]],
           'empirical_ED_v_norm_nm_inv':ed['basis_norms_nm_inv'][1],'metric_predicted_bstar_nm_inv':predicted_bstar,
           'required_ED_isotropic_scale_if_same_b':predicted_bstar/ed['basis_norms_nm_inv'][1],
           'status':'conditional relative-position compatibility; absolute ED scale and axis convention are NOT jointly validated',
           'warning':'This fits one length-scale parameter to the same three manually assigned XRD peaks. It is not a refined structure or a fit to raw data.'}
    tiff_reports=[]
    for entry in cfg['tiffs']:
        a,r=read_display_tiff(entry);tiff_reports.append(r)
        stem=Path(entry['path']).stem;write_json(out/f'{stem}_peaks.json',r)
        fig,ax=plt.subplots(figsize=(7,7));ax.imshow(a,cmap='gray',vmin=0,vmax=255)
        ps=r['peaks'];ax.scatter([p['x_px'] for p in ps],[p['y_px'] for p in ps],s=20,edgecolor='tomato',facecolors='none',lw=.6)
        cy,cx=r['center_yx_px'];ax.plot(cx,cy,'c+');ax.set_title(f"[{entry['zone_label']}] display TIFF: position candidates")
        fig.tight_layout();fig.savefig(out/f'{stem}_annotated.png',dpi=150);plt.close(fig)
    rows=[]
    for sample,key in [('A','xrd_A_resolved_d_A'),('B','xrd_B_d_A')]:
        for d in cfg[key]:rows.append({'sample':sample,'d_A':d,'reading_halfwidth_A':w,'g_nm_inv':10/d,'g_reading_interval_nm_inv':[10/(d+w),10/(d-w)],'hkl':'unassigned','method':'manual visual plot reading'})
    write_json(out/'xrd_manual_positions.json',{'source':cfg['xrd_source'],'positions':rows})
    with (out/'xrd_manual_positions.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['sample','d_A','visual_halfwidth_A','g_nm-1']);writer.writerows([[r['sample'],r['d_A'],w,r['g_nm_inv']] for r in rows])
    result={'metric_sensitivity_diagnostic_not_uncertainty':sensitivity,'conditional_relative_q_join':joint,'xrd_position_screen':xrd,'conditional_monoclinic_identity':identity,'ed_empirical_modulation':ed,
            'tiff_inputs':[{k:v for k,v in r.items() if k!='peaks'} for r in tiff_reports],
            'interpretation':'XRD positions favor quarter-b* among tested fixed-metric assignments. ED empirical fractional offset must be considered separately; do not declare a common q or a unique structure.',
            'raw_xrd':'pending','chemical_LFNSO_parent':'not supplied','atomic_structure_ranking':'not performed','gate1':'not_evaluated','gate2':'not_evaluated','stem_simulation':'not_run'}
    write_json(out/'result.json',result)
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    ds=np.linspace(0,.5,401);rms=[np.sqrt(np.mean((predicted_d(lattice,test['parent_hkl'],d)-test['B_d_A'])**2)) for d in ds]
    axes[0].plot(ds,rms,label='starting metric');
    for scale in [.99,1.01]:
        other=[np.sqrt(np.mean((scale*predicted_d(lattice,test['parent_hkl'],d)-test['B_d_A'])**2)) for d in ds]
        axes[0].plot(ds,other,ls=':',label=f'length scale {scale}')
    axes[0].legend(fontsize=8);axes[0].set(xlabel='delta along b* (conditional indexing)',ylabel='XRD position RMS error (Å)',title='Fixed starting metric / 3 manually read B features')
    vals=[abs(p['residual_uv'][1]) for p in ed['satellites']];axes[1].hist(vals,bins=np.linspace(.15,.45,31),color='slateblue')
    for ax in axes:
        for d,label in [(1/6,'1/6'),(.25,'1/4'),(1/3,'1/3')]:ax.axvline(d,color='gray',ls='--',lw=.8);ax.text(d,ax.get_ylim()[1]*.85,label,rotation=90,va='top')
    axes[1].set(xlabel='Absolute empirical residual along v (~b*)',ylabel='Detected satellite-like spots (correlated)',title='[101]-labelled DM3: empirical reciprocal basis')
    fig.tight_layout();fig.savefig(out/'position_comparison.png',dpi=160);plt.close(fig)
    notes='''<h1>LFNSO21-B: XRD 그림 + 회절 위치 비교</h1><p>그림에서 수동 판독한 XRD 위치의 범위는 ±0.02 Å로 두었습니다. 이는 통계적 표준오차가 아닙니다. 강도·피크 폭·Rietveld refinement를 수행하지 않았습니다.</p>
<p>동일한 starting metric 및 지정한 반사족을 가정하면 XRD의 세 B peak는 1/4 b*의 +1차 satellite 위치에 가깝습니다. 반면 [101] 표기 DM3의 경험적 격자에서 보이는 잔차는 별도로 평가해야 합니다. 두 결과를 강제로 한 q에 맞추지 않습니다.</p>
<p>Metric sensitivity: 시작 격자의 길이를 ±1% 바꾸는 진단에서는 최적 q가 달라집니다. ED에서 추정한 상대 q를 사용하고 길이를 약 +0.84% 보정하면 세 XRD 위치를 설명할 수 있으나, 절대 ED scale은 약 3.6%의 추가 불일치를 보여 공동 검증은 완료되지 않았습니다. 이 보정값들을 실측 calibration으로 적용하지 않았습니다.</p><p>그림의 q≈2.2, 2.48 nm⁻¹는 절대 반사 위치 1/d로 해석했습니다. 서로 다른 hkl의 반경 차이를 곧바로 3D modulation vector로 간주하지 않습니다.</p>
<p>TIFF는 8-bit 축소 display export입니다. Scale bar 기준으로 위치를 환산했으며, hkl 확정이나 약한 satellite 분리는 원본 회절·setting 확인이 필요합니다. 화학적 LFNSO 원자배열 및 구조 순위는 아직 산출하지 않았습니다.</p>'''
    body='<!doctype html><html lang="ko"><meta charset="utf-8"><title>LFNSO position screening</title><style>body{font:16px system-ui;max-width:1150px;margin:30px auto;padding:20px}pre{white-space:pre-wrap}img{max-width:100%}</style>'+notes
    for p in sorted(out.glob('*.png')):body+='<h2>'+html.escape(p.stem)+'</h2><img src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'">'
    body+='<h2>계산 결과와 가정</h2><pre>'+html.escape(json.dumps(result,indent=2,ensure_ascii=False))+'</pre></html>'
    (out/'report.html').write_text(body)
    print(json.dumps({'output':str(out),'xrd_free_delta':xrd['free_delta'],'ed_empirical_delta':ed['median_delta_v'],'ed_q_norm_nm_inv':ed['q_norm_nm_inv'],'tiff_peak_counts':[r['peak_count'] for r in tiff_reports]}))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('config');run(parser.parse_args().config)
