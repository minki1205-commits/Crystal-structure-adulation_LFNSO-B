"""Independent XRD/ED audit of the supplied C2 P/S-alternating collaborator CIF."""
import argparse,base64,csv,html,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from pymatgen.core import Lattice,Structure
from pymatgen.io.cif import CifParser
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from .physical_profile import NonnegativeProfileObjective
from .model_generator import validate
from .ed_matching import plane_reflections,register_parent,transform,in_view,match_one_to_one,electron_intensity
from .ed_metric_diagnostic import fit_correction,apply_correction
from .pipeline import digest,write_json


def predict_positions(parameters,hkl,wavelength):
    lattice=Lattice.monoclinic(*parameters[:3],parameters[3]);d=np.array([lattice.d_hkl(h) for h in hkl])
    return np.degrees(2*np.arcsin(wavelength/(2*d)))


def metric_fit(structure,observed,hkl,wavelength):
    start=np.array([structure.lattice.a,structure.lattice.b,structure.lattice.c,structure.lattice.beta])
    lower=start*np.array([.95,.95,.95,1])-np.array([0,0,0,5]);upper=start*np.array([1.05,1.05,1.05,1])+np.array([0,0,0,5])
    fit=least_squares(lambda p:predict_positions(p,hkl,wavelength)-observed,start,bounds=(lower,upper))
    trial=Structure(Lattice.monoclinic(*fit.x[:3],fit.x[3]),[site.species for site in structure],structure.frac_coords)
    return trial,{'hkl_assignments':hkl.tolist(),'observed_two_theta_deg':observed.tolist(),'predicted_two_theta_deg':predict_positions(fit.x,hkl,wavelength).tolist(),
                  'residual_two_theta_deg':fit.fun.tolist(),'RMS_position_deg':float(np.sqrt(np.mean(fit.fun**2))),
                  'parameters_a_b_c_beta':fit.x.tolist(),'length_change_percent':((fit.x[:3]/start[:3]-1)*100).tolist(),'beta_change_deg':float(fit.x[3]-start[3]),
                  'parameter_count':4,'position_count':len(observed),'nominal_residual_degrees_of_freedom':len(observed)-4,
                  'fit_bounds':[lower.tolist(),upper.tolist()],'zero_shift_deg':0.,'interpretation':'Conditional peak-position-only metric diagnostic; same fractional coordinates, no instrument/atomic refinement. 35-45 degree profile not optimized.'}


def export_calculation_copy(structure,path):
    # Reset imported labels: one disordered site needs distinct CIF rows per element.
    clean=Structure(structure.lattice,[s.species for s in structure],structure.frac_coords)
    clean.to(filename=str(path),fmt='cif',symprec=None)


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .data_io import read_dm
    from PIL import Image
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    cif=(cp.parent/cfg['cif']).resolve();original=CifParser(cif,frac_tolerance=0).parse_structures(primitive=False)[0]
    xrd_path=(cp.parent/cfg['xrd_csv']).resolve();raw=np.loadtxt(xrd_path,delimiter=',',skiprows=1);fitfile=(cp.parent/cfg['xrd_fits']).resolve();fits=json.loads(fitfile.read_text())['fits']['LFNSO21-B']
    feature_names=['anchor_001_assumed','anchor_020_assumed','feature_1','feature_2','feature_3','feature_4'];observed=np.array([fits[k]['center_two_theta_deg'] for k in feature_names])
    variants={'as_supplied':original};metric_diagnostics={}
    for name,second in [('metric_fit_110',[1,1,0]),('metric_fit_020',[0,2,0])]:
        hkl=np.array([[0,0,2],second,[0,2,1],[1,1,-2],[1,1,2],[0,2,3]])
        trial,report=metric_fit(original,observed,hkl,cfg['wavelength_A']);variants[name]=trial;metric_diagnostics[name]=report
    # Conditional identity from independent (002),(020) anchors, not a universal indexing claim.
    s002,s020=2*np.sin(np.deg2rad(observed[:2]/2));identity=np.degrees(2*np.arcsin(np.sqrt((s002/2)**2+s020**2)/2))
    summary={'config':cfg,'source_cif_sha256':digest(cif),'source_xrd_sha256':digest(xrd_path),'source_fit_sha256':digest(fitfile),'code_sha256':digest(Path(__file__)),
             'original_cell':list(original.lattice.parameters),'original_site_count':len(original),'declared_space_group':'C2','detected_space_group':SpacegroupAnalyzer(original,symprec=.01).get_space_group_symbol(),
             'composition_check':validate(original,4),'composition_note':'Li4Fe0.6667Ni0.3333SbO6 differs from .67/.33 by rounding-sized 0.0033 per formula; not treated as a structural rejection.',
             'site_hypothesis':'Ni is concentrated in the statistically mixed S layer (Li0.5 Fe0.1667 Ni0.3333), with a fully Fe P-layer site. This is distinct from uniform .67/.33 Fe/Ni co-occupancy, not automatically chemically invalid.',
             'user_zone_information':'User says collaborator STEM matching used [001] and [100] of this CIF. Existing ED file labels remain [001], [101], [110]. Their correspondence is tested as hypotheses.',
             'conditional_identity':{'anchors_002_and_020_deg':observed[:2].tolist(),'predicted_021_deg':float(identity),'observed_feature1_deg':float(observed[2]),'residual_deg':float(identity-observed[2]),'note':'Only if the ~19.834 degree feature is (020). It can instead be tested as the strong (110); that alternate assignment is explicitly fit separately.'},
             'metric_diagnostics':metric_diagnostics,'xrd':{},'ed':[],'validated_structure':False}
    windows={'superstructure':(20.4,34.6),'high_angle_postfit':(35,45),'low_angle':(17,20.3),'full_10_80':(10,80)}
    patterns={};profile_figures=[]
    for name,structure in variants.items():
        folder=out/name;folder.mkdir();export_calculation_copy(structure,folder/'calculation_copy.cif')
        p=XRDCalculator(wavelength=cfg['wavelength_A']).get_pattern(structure,two_theta_range=(9,81));patterns[name]=p
        record={'two_theta_deg':p.x.tolist(),'relative_intensity':p.y.tolist(),'hkls':p.hkls,'d_A':np.asarray(p.d_hkls).tolist()};write_json(folder/'xrd_reflections.json',record)
        scores={}
        for width in [.15,.25,.4]:
            scores[str(width)]={}
            for win,(lo,hi) in windows.items():
                mask=(raw[:,0]>=lo)&(raw[:,0]<=hi);r=NonnegativeProfileObjective(p.x,*raw[mask].T,width).evaluate(p.y)
                scores[str(width)][win]={k:v for k,v in r.items() if k!='profile'}
        summary['xrd'][name]={'metric':list(structure.lattice.parameters),'scores':scores,'geometry':{k:v for k,v in validate(structure,4).items() if k!='coordination'}}
        fig,axes=plt.subplots(1,3,figsize=(15,4))
        for ax,(win,(lo,hi)) in zip(axes,list(windows.items())[:3]):
            mask=(raw[:,0]>=lo)&(raw[:,0]<=hi);x,y=raw[mask].T;r=NonnegativeProfileObjective(p.x,x,y,.25).evaluate(p.y)
            ax.plot(x,y,'k-',lw=.8,label='LFNSO21-B measured');ax.plot(x,r['profile'],'r-',lw=.8,label='CIF prediction + nonnegative baseline');ax.set(xlabel='2theta (degrees)',title=f"{win}: R2={r['R2_vs_nonnegative_linear_background']:.3f}");ax.legend(fontsize=7)
        fig.suptitle(name+' / Cu Kalpha 1.5406 A provisional; FWHM 0.25 degrees');fig.tight_layout();f=name+'_xrd.png';fig.savefig(out/f,dpi=140);plt.close(fig);profile_figures.append(f)
    original_hkl=np.array([[0,0,2],[1,1,0],[0,2,1],[1,1,-2],[1,1,2],[0,2,3]])
    original_params=np.array([original.lattice.a,original.lattice.b,original.lattice.c,original.lattice.beta]);initial_positions=predict_positions(original_params,original_hkl,cfg['wavelength_A'])
    summary['as_supplied_six_feature_comparison']={'hkl_provisional':original_hkl.tolist(),'observed_two_theta_deg':observed.tolist(),'predicted_two_theta_deg':initial_positions.tolist(),'predicted_minus_observed_deg':(initial_positions-observed).tolist()}
    edp=(cp.parent/cfg['ed_summary']).resolve();ed=json.loads(edp.read_text());summary['source_ed_summary_sha256']=digest(edp)
    zone_sets={'DP01_101':[[1,0,1],[1,0,-1],[5,0,3],[7,0,-3],[0,0,1],[1,0,0]],'dp3_110':[[1,1,0],[1,0,0]],'OneView_300X_200kV_0021_001':[[0,0,1],[-1,0,3]],'OneView_250X_200kV_0052':[[1,1,0],[1,0,0]]}
    summary['mapped_zone_note']='Additional mapped hypotheses use c_2layer = 2*c_reference + a/3, inferred from idealized metrics; not a confirmed experimental orientation or integer basis equivalence. [001]->[-103], [101]->[503], [10-1]->[70-3]. Literal CIF zones remain primary where supplied by user.'
    overlays=[]
    for name in ['as_supplied','metric_fit_110']:
        structure=variants[name]
        for dataset in ed['datasets']:
            if dataset['id'] in cfg.get('raw_overrides',{}):
                dataset['raw']=str((cp.parent/cfg['raw_overrides'][dataset['id']]).resolve())
            if digest(Path(dataset['raw']))!=dataset['raw_sha256']:raise ValueError('ED image hash mismatch: '+dataset['id'])
            for zone in zone_sets[dataset['id']]:
                key=name+'__'+dataset['id']+'__'+','.join(map(str,zone));folder=out/name
                hk,hp,g,xy=plane_reflections(structure.lattice,1,zone,14);reference=xy[(hk[:,0]+hk[:,1])%2==0]
                reg=register_parent(np.array(dataset['observed']),reference,[.9,1.1],3*dataset['pixel_scale'])
                reference_image=transform(reference,reg['parameters'],reg['mirror'])
                try:correction=fit_correction(dataset['observed'],reference_image,reg['parameters'][2:],3*dataset['pixel_scale'],.35 if dataset['kind']=='dm3' else .48,.05 if dataset['kind']=='dm3' else .16)
                except ValueError as e:correction={'status':'not_fitted','reason':str(e)}
                pred=transform(xy,reg['parameters'],reg['mirror']);intensity=electron_intensity(structure,hk,np.linalg.norm(g,axis=1));scores=[]
                write_json(folder/(dataset['id']+'__'+','.join(map(str,zone))+'_reflections.json'),{'hkl_CIF':hk.tolist(),'g_similarity_nm_inv':pred.tolist(),'kinematic_intensity':intensity.tolist()})
                for mode in ['similarity','affine_diagnostic']:
                    if mode=='affine_diagnostic' and correction.get('status')=='not_fitted':continue
                    points=pred if mode=='similarity' else apply_correction(pred,correction);visible=in_view(points,dataset);norm=max(float(intensity[visible].max()) if visible.any() else 0,1e-30)
                    for threshold in [.1,1,5]:
                        selected=visible&(intensity>=norm*threshold/100);ids=np.flatnonzero(selected)
                        for tolerance in dataset['position_tolerances_nm_inv']:
                            m=match_one_to_one(dataset['observed'],points[selected],tolerance)
                            for pair in m['pairs']:pair['source_peak_index']=dataset['source_peak_indices'][pair['observed_index']];pair['reflection_index']=int(ids[pair['predicted_index']])
                            scores.append({'mode':mode,'threshold_percent':threshold,'tolerance_nm_inv':tolerance,**m})
                grid_scores=[]
                for mode in ['similarity','affine_diagnostic']:
                    if mode=='affine_diagnostic' and correction.get('status')=='not_fitted':continue
                    points=pred if mode=='similarity' else apply_correction(pred,correction)
                    for label,allowed in [('integer_grid',np.ones(len(hk),bool)),('C_centered_grid',(hk[:,0]+hk[:,1])%2==0)]:
                        selected=in_view(points,dataset)&allowed
                        m=match_one_to_one(dataset['observed'],points[selected],dataset['position_tolerances_nm_inv'][1])
                        grid_scores.append({'mode':mode,'grid':label,**{k:v for k,v in m.items() if k!='pairs'}})
                row={'variant':name,'dataset':dataset['id'],'nominal_ED_zone':dataset['nominal_zone'],'CIF_zone_hypothesis':zone,'registration':reg,'affine_diagnostic':correction,'scores':scores,'intensity_free_grid_diagnostic':grid_scores};summary['ed'].append(row)
                standard=[r for r in scores if r['threshold_percent']==1 and r['tolerance_nm_inv']==dataset['position_tolerances_nm_inv'][1]]
                print(key,'scale',round(reg['parameters'][1],4),'parent RMS',round(reg['trimmed_RMS_nm_inv'],4),'F1',[(r['mode'],round(r['f1'],3)) for r in standard],flush=True)
                # Only supplied CIF literal primary/projection alternatives are drawn, not a best-zone cherry pick.
                if name=='as_supplied' and zone in [[0,0,1],[1,1,0],[1,0,0]]:
                    if dataset['kind']=='dm3':a,_=read_dm(Path(dataset['raw']));a=np.log1p(np.maximum(a,0));lo,hi=np.percentile(a,[1,99.8])
                    else:a=np.array(Image.open(dataset['raw']).convert('L'));lo,hi=0,255
                    mode='affine_diagnostic' if correction.get('status')!='not_fitted' else 'similarity';pp=apply_correction(pred,correction) if mode=='affine_diagnostic' else pred
                    visible=in_view(pp,dataset);selected=visible&(intensity>=max(float(intensity[visible].max()),1e-30)*.01);pixels=pp[selected]/dataset['pixel_scale']+dataset['center'];obs=np.array(dataset['observed'])/dataset['pixel_scale']+dataset['center']
                    fig,ax=plt.subplots(figsize=(7,7));ax.imshow(a,cmap='gray',vmin=lo,vmax=hi);ax.scatter(obs[:,1],obs[:,0],s=40,edgecolors='tomato',facecolors='none',lw=.7,label='Measured detections');ax.scatter(pixels[:,1],pixels[:,0],s=20,marker='+',c='cyan',lw=.7,label='Supplied CIF >=1% kinematic')
                    y0,x0,y1,x1=dataset['roi'];ax.set(xlim=(x0,x1),ylim=(y1,y0),title=f"{dataset['id']} vs CIF {zone}\n{mode}; axis hypothesis, not reassigned file label");ax.legend(fontsize=7);fig.tight_layout();f=key+'.png';fig.savefig(out/f,dpi=110);plt.close(fig);overlays.append(f)
    summary['limitations']=['CIF comments are authors assumptions, not measurements or instructions.','The two metric fits test alternate reflection assignments; they are not Rietveld refinements or validations.','Experimental STEM images exist, but the collaborator simulated images, image-to-zone registration, thickness and detector parameters were not supplied; their claimed match was not independently validated.','Electron factors include all site occupancies. ED intensity thresholds are diagnostics; no dynamical scattering, tilt, thickness, or measured intensity fit.','New CIF-specific pose is needed because its cell differs. All variants get the same registration algorithm/ranges, but not a validated common calibration.','DP01 registration itself has large residuals; zero matches in tested settings is not proof against all possible orientations.','Measured detections may contain artifacts; all unmatched detections are not automatically physical satellites.','Partial agreement in one ED view cannot establish the complete XRD/ED/STEM structure.']
    write_json(out/'summary.json',summary)
    with (out/'ed_scores.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['variant','dataset','CIF_zone','mode','threshold_percent','tolerance_nm_inv','matches','observed','predicted','recall','precision','F1','RMS_nm_inv'])
        for r in summary['ed']:
            for v in r['scores']:w.writerow([r['variant'],r['dataset'],str(r['CIF_zone_hypothesis']),v['mode'],v['threshold_percent'],v['tolerance_nm_inv'],v['matches'],v['observed_count'],v['predicted_count'],v['recall'],v['precision'],v['f1'],v['matched_RMS_nm_inv']])
    page='<!doctype html><meta charset="utf-8"><title>Collaborator C2 P/S CIF audit</title><style>body{font:16px system-ui;max-width:1300px;margin:30px auto}img{max-width:100%}pre{white-space:pre-wrap}td,th{padding:7px;border-bottom:1px solid #ccc}</style><h1>협력기관 C2 P/S 교대 CIF · 독립 XRD/ED 비교</h1>'
    page+='<p>원본을 보존하고 받은 셀 그대로 계산한 결과와, 두 가지 hkl 가정으로 metric만 조건부 조정한 진단을 구분했습니다. 사용자에 따르면 협력기관 STEM은 이번 CIF의 [001]/[100] 기준입니다. 기존 [110] ED와 CIF [100]을 같은 방향으로 확정하지 않고 별도 가설로 시험했습니다.</p>'
    page+='<table><tr><th>조건</th><th>superstructure R²</th><th>35–45° R²</th><th>17–20.3° R²</th></tr>'
    for name,r in summary['xrd'].items():page+='<tr><td>'+name+'</td>'+''.join(f'<td>{r["scores"]["0.25"][w]["R2_vs_nonnegative_linear_background"]:.4f}</td>' for w in ['superstructure','high_angle_postfit','low_angle'])+'</tr>'
    page+='</table><p>R²는 비음수 선형 background-only 대비 개선량이며 정확도나 구조 확률이 아닙니다. Metric fit은 여섯 저각 feature와 네 격자 변수를 사용하며 35–45°는 목적함수에 넣지 않았습니다. 피크 폭 .15/.25/.40°도 따로 평가했습니다.</p>'
    page+='<p><b>판정: 받은 CIF 그대로는 XRD 전체를 설명하지 못합니다. 다만 [110] ED에는 부분적으로 좋은 대응이 있어 모델 전체를 즉시 배제할 근거로 삼지 않습니다.</b> 저각 피크만 맞춘 격자는 35–45° 구간에서 악화됩니다. [001]·DP01까지 함께 설명하는지는 미해결이며 STEM 정합 주장은 별도 검증이 필요합니다.</p>'
    page+='<table><tr><th>원본 CIF / 영상</th><th>CIF zone</th><th>기본 위치 대응</th><th>추가 affine 대응</th><th>관측 점</th><th>기본 배율</th></tr>'
    for row in summary['ed']:
        if row['variant']!='as_supplied':continue
        a=[v for v in row['scores'] if v['threshold_percent']==1 and v['mode']=='similarity'][1]
        b=[v for v in row['scores'] if v['threshold_percent']==1 and v['mode']=='affine_diagnostic']
        page+=f'<tr><td>{row["dataset"]}</td><td>{row["CIF_zone_hypothesis"]}</td><td>{a["matches"]}</td><td>{b[1]["matches"] if b else "N/A"}</td><td>{a["observed_count"]}</td><td>{row["registration"]["parameters"][1]:.4f}</td></tr>'
    page+='</table><p>강도 문턱 1%, 위치 허용오차 DM3 0.08 / TIFF 0.24 nm⁻¹. 배율은 각 영상에서 ±10% 허용한 위치 정합값이며 검증된 calibration 보정이 아닙니다. DP01은 정합 자체의 잔차가 크므로 시험하지 않은 축까지 배제할 수 없습니다. TIFF는 표시용 영상이며 검출 점에 artifact가 포함될 수 있습니다.</p>'
    page+='<table><tr><th>원본 CIF / 영상 / zone</th><th>0.1% 대응 수</th><th>1% 대응 수</th><th>5% 대응 수</th><th>C격자 위치만 / 관측 수</th></tr>'
    for row in summary['ed']:
        if row['variant']!='as_supplied' or row['affine_diagnostic'].get('status')=='not_fitted':continue
        counts=[[v for v in row['scores'] if v['mode']=='affine_diagnostic' and v['threshold_percent']==t][1]['matches'] for t in [.1,1,5]]
        grid=next(v for v in row['intensity_free_grid_diagnostic'] if v['mode']=='affine_diagnostic' and v['grid']=='C_centered_grid')
        page+=f'<tr><td>{row["dataset"]} / {row["CIF_zone_hypothesis"]}</td>'+''.join(f'<td>{n}</td>' for n in counts)+f'<td>{grid["matches"]} / {grid["observed_count"]}</td></tr>'
    page+='</table><p>위 표는 affine 진단 기준입니다. [001]은 강도 문턱에 민감하며 격자 위치 자체의 대응은 더 높습니다. 따라서 1% 기준의 미대응 점을 모두 격자 불일치나 구조 반증으로 볼 수 없습니다.</p>'
    for f in profile_figures+overlays:page+='<img src="data:image/png;base64,'+base64.b64encode((out/f).read_bytes()).decode()+'">'
    page+='<pre>'+html.escape(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre>'
    (out/'report.html').write_text(page);(out/'source_snapshot.py').write_bytes(Path(__file__).read_bytes());return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
