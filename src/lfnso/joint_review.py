"""Summarize XRD/ED tensions without declaring a joint structure winner."""
import argparse,base64,csv,html,json
from pathlib import Path
import numpy as np
from .pipeline import digest,write_json


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    paths={key:(cp.parent/cfg[key]).resolve() for key in ['ed','affine','positions','xrd']};sources={key:json.loads(p.read_text()) for key,p in paths.items()}
    ed=sources['ed'];affine=sources['affine'];q=sources['positions']['ed_empirical_modulation'];xrd=sources['xrd']['conditional_pure_bstar_modulation']
    values=np.array([abs(p['residual_uv'][1]) for p in q['satellites']]);qrows=[]
    for period,harmonic in [(4,1),(7,2),(11,3)]:
        delta=harmonic/period;err=values-delta
        qrows.append({'period_b':period,'harmonic':harmonic,'delta':delta,'ED_offset_RMS_fraction_v':float(np.sqrt(np.mean(err**2))),
                      'ED_median_absolute_offset_fraction_v':float(np.median(abs(err))),'XRD_conditional_delta_minus_harmonic':xrd['delta']-delta})
    dataset_order=['DP01_101','OneView_300X_200kV_0021_001','dp3_110','OneView_250X_200kV_0052']
    matrix=[]
    for model in ed['models']:
        if model['role']!='trial':continue
        row={'model':model['id'],'metric':model['metric'],'period_b':model['period'],'XRD_R2_20p4_34p6':model['xrd_R2_full_superstructure'],'ED_F1_similarity':{},'ED_F1_affine_diagnostic':{}}
        for dataset in dataset_order:
            e=next(r for r in ed['rows'] if r['model']==model['id'] and r['dataset']==dataset and (dataset!='DP01_101' or r['zone_hypothesis']==[1,0,-1]))
            a=next(r for r in affine['rows'] if r['pattern_key']==e['pattern_key'])
            row['ED_F1_similarity'][dataset]=e['standard_score']['f1'];row['ED_F1_affine_diagnostic'][dataset]=a['affine_standard_score']['f1']
        matrix.append(row)
    summary={'config':cfg,'source_sha256':{k:digest(p) for k,p in paths.items()},'code_sha256':digest(Path(__file__)),
             'process':{'data_inventory':'completed for accessible uploads; microscope metadata incomplete','hypothesis_generation':'32 earlier library models plus 20 occupancy-search endpoints; unrelaxed','XRD_screening':'completed exploratory comparisons; no gate pass','ED_position_matching':'12 latest candidates + 2 ordered controls, 4 images, 70 model/image/zone comparisons; similarity plus bounded affine diagnostic','joint_structure_selection':'unresolved','STEM_multislice':'not run'},
             'q_comparison':{'ED_median_fraction_v':q['median_delta_v'],'ED_correlated_satellite_count':len(values),'XRD_conditional_fraction_bstar':xrd['delta'],'hypotheses':qrows,'note':'ED values reuse the same DP01 data; no extra independent evidence. Conditional indexing and v~b*; no commensurability or common 3D q is established.'},
             'score_vectors':matrix,'validated_shortlist':[],
             'conclusions':['The 11b XRD-best trial is not a joint XRD/ED winner.',
                            'Conditional ED fractional spacing is closer to 2/7 than 3/11, but candidate harmonics/intensities and all zones still disagree.',
                            '[10-1] fits the [101]-named image better within the tested setting. Filename and user zone label are preserved.',
                            'Calibration/metric distortions materially affect matching, and the [001] similarity scale reaches the tested lower bound for the starting metric.',
                            'In the exact [110] ZOLZ, b-only supercell size does not produce new fractional positions; occupancy changes intensities only there.',
                            'Next: establish indexing/calibration and explore limited coupled c/stacking/displacement hypotheses; do not optimize a single ED overlay at the expense of XRD.']}
    write_json(out/'summary.json',summary)
    with (out/'score_vectors.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['model','period_b','metric','XRD_R2']+['similarity_'+d for d in dataset_order]+['affine_diagnostic_'+d for d in dataset_order])
        for r in matrix:w.writerow([r['model'],r['period_b'],r['metric'],r['XRD_R2_20p4_34p6']]+list(r['ED_F1_similarity'].values())+list(r['ED_F1_affine_diagnostic'].values()))
    fig,ax=plt.subplots(figsize=(9,4))
    ax.scatter(values,np.arange(len(values)),color='black',s=20,label='Selected DP01 satellite offsets (correlated)')
    for r,color in zip(qrows,['gray','tab:blue','tab:orange']):ax.axvline(r['delta'],color=color,label=f"{r['harmonic']}/{r['period_b']} b*")
    ax.axvline(xrd['delta'],color='tab:red',ls='--',label='Conditional XRD delta')
    ax.set(xlabel='Fractional offset: empirical v for ED; assumed b* for XRD',ylabel='Selected spot index',title='Relative spacing diagnostic; indexing and scale remain conditional');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'q_comparison.png',dpi=160);plt.close(fig)
    page='<!doctype html><meta charset="utf-8"><title>LFNSO21-B XRD + ED progress</title><style>body{font:16px system-ui;max-width:1250px;margin:30px auto}img{max-width:100%}td,th{padding:7px;border-bottom:1px solid #ccc}pre{white-space:pre-wrap}</style><h1>LFNSO21-B · XRD와 ED 공동 검토</h1>'
    page+='<p><b>현재: 후보 생성과 XRD screening 이후, 실제 ED 위치 정합까지 진행했습니다. 최종 구조 선택은 미완료입니다.</b> 최신 후보 12개와 ordered control 2개를 DM3 두 장·TIFF 두 장에 대조했습니다. [101]과 [10-1] 대안을 포함해 70개 조합, 세 반사 threshold·세 위치 허용오차를 평가했습니다. 추가 방향별 보정은 원인 분리를 위한 민감도 진단입니다.</p>'
    page+='<p>지금까지 ED는 경험적 q 및 기하학적 표 중심이었고, 이번에 실제 검출 좌표와 후보 패턴의 일대일 매칭을 추가했습니다. Fe/Ni 점유율을 반영한 전자 산란인자를 사용하지만 두께·tilt·다중 산란은 계산하지 않았습니다. STEM은 아직 실행하지 않았습니다.</p>'
    page+='<p>조건부 XRD 최상 11b 모델이 ED까지 동시에 최상인 것은 아닙니다. ED 위성의 상대 간격은 2/7 쪽에 가깝고, 단순 b방향 변조만으로는 정확한 [110] 평면에 새로운 fractional 반사를 만들 수 없습니다. 이것이 stacking/c축 결합 또는 변위 가설을 분리해서 시험해야 하는 이유입니다. 현재 데이터만으로 그 원인을 확정하지 않습니다.</p>'
    page+='<img src="data:image/png;base64,'+base64.b64encode((out/'q_comparison.png').read_bytes()).decode()+'">'
    page+='<h2>후보별 점수 벡터 · 합산 순위 아님</h2><p>R²는 XRD 선형 배경 대비 개선량, F1은 ED 위치 매칭의 precision/recall 조화평균입니다. 서로 다른 정의의 점수를 합산하지 않았습니다. 아래 ED는 회전+등방 배율 모드, threshold 1%, DM3 0.08 nm⁻¹ / TIFF 0.24 nm⁻¹입니다. DP01 열은 미확정 [10-1] 가설입니다.</p><table><tr><th>후보</th><th>XRD R²</th><th>DP01 F1</th><th>[001] TIFF F1</th><th>[110] DM3 F1</th><th>[110] TIFF F1</th></tr>'
    for row in sorted(matrix,key=lambda r:r['XRD_R2_20p4_34p6'],reverse=True):page+='<tr><td>'+row['model']+'</td>'+''.join(f'<td>{v:.3f}</td>' for v in [row['XRD_R2_20p4_34p6'],*row['ED_F1_similarity'].values()])+'</tr>'
    page+='</table><h2>검증 범위와 다음 단계</h2><p>원본·CIF는 보존했습니다. 큰 reciprocal grid의 우연 매칭을 확인하기 위해 동일 모델 회전 대조군과 precision을 함께 기록했고, 반사 threshold/위치 허용오차 민감도를 저장했습니다. 미검출 반사를 모두 실제 소멸로 간주하지 않습니다. 결정축 convention과 영상 calibration 확인, XRD metric/profile 검토, c축·적층·변위의 제한 가설 확장이 남아 있습니다.</p><pre>'+html.escape(json.dumps(summary,indent=2,ensure_ascii=False))+'</pre>'
    (out/'report.html').write_text(page);(out/'joint_review_snapshot.py').write_bytes(Path(__file__).read_bytes());return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
