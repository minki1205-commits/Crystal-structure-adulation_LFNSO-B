"""Review c-layer controls and direct joint occupational searches without a gate claim."""
import argparse,base64,csv,html,json
from pathlib import Path
import numpy as np
from pymatgen.io.cif import CifParser
from .occupational_search import ProfileObjective
from .coupled_search import powder,verify_powder
from .pipeline import digest,write_json
from .physical_profile import NonnegativeProfileObjective


def edscore(row,dataset,mode='affine_diagnostic'):
    return next(r['f1'] for r in row['ed_scores'] if r['dataset']==dataset and r['mode']==mode and r['threshold_percent']==1)


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    staged_path=(cp.parent/cfg['staged_summary']).resolve();staged=json.loads(staged_path.read_text());rows=[];locations={};source_hashes={str(staged_path):digest(staged_path)}
    for item in cfg['joint_summaries']:
        path=(cp.parent/item).resolve();source_hashes[str(path)]=digest(path)
        for row in json.loads(path.read_text())['rows']:rows.append(row);locations[row['id']]=path.parent/row['id']
    legacy_records=[]
    for item in cfg.get('legacy_joint_summaries',[]):
        path=(cp.parent/item).resolve();source_hashes[str(path)]=digest(path)
        legacy=json.loads(path.read_text())
        raw_legacy=np.loadtxt(cfg['xrd_csv'],delimiter=',',skiprows=1);mask=(raw_legacy[:,0]>=20.4)&(raw_legacy[:,0]<=34.6)
        for row in legacy['rows']:
            pat=json.loads((path.parent/row['id']/'xrd.json').read_text())
            v=NonnegativeProfileObjective(np.array(pat['two_theta_deg']),*raw_legacy[mask].T,.25).evaluate(np.array(pat['intensity_arbitrary_common_units']))
            legacy_records.append({'cycle':path.parent.name,'id':row['id'],'old_unconstrained_R2':row['xrd_scores']['superstructure']['R2_vs_linear_background'],'same_CIF_nonnegative_background_R2':v['R2_vs_nonnegative_linear_background'],'old_linear_background':row['xrd_scores']['superstructure']['linear_background']})
    best_xrd=max(rows,key=lambda r:r['xrd_scores']['superstructure']['R2_vs_linear_background'])
    best_110=max(rows,key=lambda r:(edscore(r,'dp3_110'),r['xrd_scores']['superstructure']['R2_vs_linear_background']))
    compromise=next(r for r in rows if r['id']==cfg['illustrative_compromise'])
    selected=[best_xrd,compromise,best_110]
    # Independent numerical verification on the exported selected CIFs.
    verification={}
    for row in selected:
        path=locations[row['id']]/'model.cif';s=CifParser(path,frac_tolerance=0).parse_structures(primitive=False)[0]
        aa,ii=powder(s);v=verify_powder(s,aa,ii,1.5406)
        saved=json.loads((locations[row['id']]/'xrd.json').read_text());np.testing.assert_allclose(aa,saved['two_theta_deg'],atol=1e-7)
        previous=np.array(saved['intensity_arbitrary_common_units']);np.testing.assert_allclose(ii/ii.max(),previous/previous.max(),atol=3e-6)
        verification[row['id']]={'cif_sha256':digest(path),'CIF_roundtrip_powder_verified':True,**v}
    summary={'config':cfg,'source_sha256':source_hashes,'review_code_sha256':digest(Path(__file__)),'profile_code_sha256':digest(Path(__file__).with_name('physical_profile.py')),
             'staged_trial_count':staged['candidate_count'],'staged_geometry_eligible':staged['geometry_eligible_count'],'direct_search_endpoints':len(rows),
             'direct_search_swap_proposals':sum(r['proposal_steps'] for r in rows),'direct_search_nontrivial_evaluations':sum(r['nontrivial_evaluations'] for r in rows),
             'legacy_unconstrained_background_audit':legacy_records,'selected_for_illustration_not_shortlist':[r['id'] for r in selected],'independent_numerical_verification':verification,
             'rows':rows,'validated_shortlist':[],
             'conclusions':['Layer sequence/registry trials alone did not solve the joint pattern.','Direct six-c within-layer permutations substantially improved exploratory XRD fit.','Increasing [110] ED support weight traded XRD fit for additional position matches, but the [101]-named image under the alternate [10-1] setting worsened.','No structure is established; six-c periodicity, oxidation states and dynamical diffraction remain unresolved.','The next discriminating constraint is maintaining in-plane satellite ordering across [001]/[101] while explaining [110], with calibrated zone/metric handling.']}
    write_json(out/'summary.json',summary)
    raw=np.loadtxt(cfg['xrd_csv'],delimiter=',',skiprows=1);base_pattern=json.loads((cp.parent/cfg['baseline_xrd']).resolve().read_text())
    figures=[]
    for row in selected:
        pattern=json.loads((locations[row['id']]/'xrd.json').read_text());fig,axes=plt.subplots(1,3,figsize=(15,4))
        for ax,(name,lo,hi) in zip(axes,[('superstructure',20.4,34.6),('high_angle',35,45),('low_angle',17,20.3)]):
            mask=(raw[:,0]>=lo)&(raw[:,0]<=hi);x,y=raw[mask].T
            fit=NonnegativeProfileObjective(np.array(pattern['two_theta_deg']),x,y,.25).evaluate(np.array(pattern['intensity_arbitrary_common_units']))
            baseline=NonnegativeProfileObjective(np.array(base_pattern['two_theta_deg']),x,y,.25).evaluate(np.array(base_pattern['intensity_arbitrary_common_units']))
            ax.plot(x,y,'k-',lw=.7,label='B measured');ax.plot(x,baseline['profile'],color='steelblue',lw=.7,label='Earlier 11b trial');ax.plot(x,fit['profile'],'r-',lw=.8,label='New 3D trial')
            ax.set(title=f"{name}: R2={fit['R2_vs_nonnegative_linear_background']:.3f}",xlabel='2theta (degrees)');ax.legend(fontsize=7)
        fig.suptitle(row['id']+'; fitted data, not independent validation');fig.tight_layout();name=row['id']+'_xrd.png';fig.savefig(out/name,dpi=140);plt.close(fig);figures.append(name)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for row in rows:
        period=row['base']['b_repeat'];color='tab:blue' if period==7 else 'tab:orange';x=row['xrd_scores']['superstructure']['R2_vs_linear_background']
        for ax,d in zip(axes,['dp3_110','DP01_101']):
            ax.scatter(x,edscore(row,d),c=color,s=30+20*row['ed_objective_weight']);ax.annotate(str(row['ed_objective_weight']),(x,edscore(row,d)),fontsize=7,xytext=(3,3),textcoords='offset points')
    for ax,d in zip(axes,['[110] DM3 (optimized support)','DP01 [10-1] hypothesis (post-search check)']):ax.set(xlabel='XRD superstructure R2 vs linear background',ylabel='ED position F1: affine diagnostic / 1% threshold',title=d)
    fig.suptitle('Blue: 7b x 6c; orange: 11b x 6c; labels: ED objective weight');fig.tight_layout();fig.savefig(out/'tradeoffs.png',dpi=150);plt.close(fig)
    with (out/'scores.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['model','ED_weight','seed','XRD_R2_super','XRD_R2_high','ED110_F1_affine','DP01_F1_affine','ED001_F1_affine','ED110_TIFF_F1_affine','ED110_F1_similarity'])
        for r in rows:w.writerow([r['id'],r['ed_objective_weight'],r['seed'],r['xrd_scores']['superstructure']['R2_vs_linear_background'],r['xrd_scores']['high_angle']['R2_vs_linear_background'],edscore(r,'dp3_110'),edscore(r,'DP01_101'),edscore(r,'OneView_300X_200kV_0021_001'),edscore(r,'OneView_250X_200kV_0052'),edscore(r,'dp3_110','similarity')])
    page='<!doctype html><meta charset="utf-8"><title>LFNSO21-B 3D joint search</title><style>body{font:16px system-ui;max-width:1300px;margin:30px auto}img{max-width:100%}td,th{padding:7px;border-bottom:1px solid #ccc}pre{white-space:pre-wrap}</style><h1>3D 가설 확장 및 XRD–ED 직접 탐색</h1>'
    page+=f'<p><b>아직 최종 구조를 찾지 못했습니다. 그러나 6c 층별 배열을 직접 바꾸면 XRD 피팅이 개선되는 후보를 얻었습니다.</b> 정해진 층 순서·registry·변위 시험 {staged["candidate_count"]}개(기존/반복 대조 포함)와 직접 탐색 종점 {len(rows)}개를 추가했습니다. 최종 비음수 배경 탐색은 {summary["direct_search_swap_proposals"]:,} proposal, {summary["direct_search_nontrivial_evaluations"]:,}회 비동일 원소 교환 평가를 수행했습니다. 서로 다른 구조 개수 또는 완전 탐색을 뜻하지 않습니다.</p>'
    page+='<p>Li4Fe0.67Ni0.33SbO6, Fe/Ni 공동 점유와 층별 Li:Sb:Fe-host 개수를 보존했습니다. 산소·interlayer Li 좌표와 metric은 고정했습니다. 큰 b/6 registry 시험 세 개는 CN6 비율 0.75로 제외했습니다. 나머지 기하 검사는 에너지 안정성 또는 정확한 octahedral bond-angle 검증이 아닙니다.</p>'
    page+='<p><b>배경 제약 재검토:</b> 앞선 16개 진단 종점에서 최고 0.849는 음수 선형 배경을 허용한 값이었습니다. 같은 CIF에 비음수 배경을 강제하면 약 0.783입니다. 이를 최종 개선으로 사용하지 않고, 배경 양끝과 구조 scale을 모두 비음수로 제한하여 새로 탐색한 16개 종점을 아래에 제시합니다. 새 후보의 점수는 다른 원자 배열에서 얻은 값입니다.</p>'
    page+='<p>비음수 배경 XRD R²의 제곱근과 [110] DM3의 기하학적 위치 support를 직접 최적화했습니다. ED 가중치 0/0.4/1.2/3.0은 진단 선택이며 통계적 likelihood가 아닙니다. 실제 electron intensity를 피팅한 것이 아닙니다. 기존 affine 진단 pose를 그대로 사용했고 다른 ED 방향은 사후 평가했습니다. 이미 살펴본 자료이므로 독립 검증으로 부르지 않습니다.</p>'
    page+='<table><tr><th>설명용 후보</th><th>XRD super R²</th><th>[110] DM3 F1</th><th>DP01 F1 ([10-1] 가설)</th></tr>'
    page+='<tr><td>이전 conditional B 11b</td><td>0.705</td><td>0.239</td><td>0.182</td></tr>'
    for title,r in zip(['새 XRD 최상','XRD/ED tradeoff 예시','새 [110] F1 최상'],selected):page+=f'<tr><td>{title}: {r["id"]}</td><td>{r["xrd_scores"]["superstructure"]["R2_vs_linear_background"]:.3f}</td><td>{edscore(r,"dp3_110"):.3f}</td><td>{edscore(r,"DP01_101"):.3f}</td></tr>'
    page+='</table><p>R²는 비음수 선형 background-only fit 대비 개선량이며 Rwp가 아닙니다. F1은 1% simulated-intensity threshold와 DM3 0.08 nm⁻¹ 위치 허용오차의 매칭 점수입니다. 구조 확률로 읽거나 두 값을 합쳐 최종 순위로 삼지 않습니다. 가장 잘 맞은 후보도 다른 ED 방향의 불일치가 남습니다.</p>'
    for name in ['tradeoffs.png']+figures:page+='<img src="data:image/png;base64,'+base64.b64encode((out/name).read_bytes()).decode()+'">'
    page+='<p>다음 핵심은 [001]/[101]의 in-plane satellite ordering을 유지하면서 [110] 반사를 함께 설명하는 것입니다. Zone convention·calibration과 oxygen/metric/displacement의 영향도 미해결입니다. 6c, 특정 적층 순서, 산화수, 산소 결손 또는 STEM 검증 통과를 확정하지 않았습니다.</p><pre>'+html.escape(json.dumps(summary,indent=2,ensure_ascii=False))+'</pre>'
    (out/'report.html').write_text(page);(out/'coupled_review_snapshot.py').write_bytes(Path(__file__).read_bytes());return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
