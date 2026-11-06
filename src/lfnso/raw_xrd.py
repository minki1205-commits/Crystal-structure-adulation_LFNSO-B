"""Read measured profiles and fit local peak envelopes without inventing a wavelength."""
import argparse,base64,csv,html,json
from pathlib import Path
import numpy as np
from openpyxl import load_workbook
from scipy.optimize import least_squares
from scipy.signal import find_peaks,savgol_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .pipeline import digest,write_json
from .xrd_hypotheses import compare
from .peak_groups import fit_group,group_profile


def read_workbook(path,sheet):
    wb=load_workbook(path,read_only=True,data_only=True)
    try:
        ws=wb[sheet];rows=list(ws.iter_rows(values_only=True));samples={};metadata=[]
        for column in range(0,len(rows[0]),2):
            name=rows[0][column]
            if name is None:continue
            if '2' not in str(rows[1][column]) or 'Intensity' not in str(rows[1][column+1]):
                raise ValueError('Expected angle/intensity headers')
            pairs=[]
            for row_number,row in enumerate(rows[2:],3):
                pair=row[column:column+2]
                if len(pair)!=2 or any(v is None for v in pair):raise ValueError(f'Missing data at Excel row {row_number}, sample {name}')
                if not all(isinstance(v,(int,float)) for v in pair):raise ValueError('Non-numeric or uncached formula in workbook')
                pairs.append(pair)
            a=np.array(pairs,float)
            if not np.isfinite(a).all() or (np.diff(a[:,0])<=0).any():raise ValueError('Finite, strictly increasing angle data required')
            samples[name]=a
            metadata.append({'sample':name,'excel_columns_1_based':[column+1,column+2], 'rows_1_based':[3,len(rows)],
                             'points':len(a),'two_theta_range_deg':[float(a[0,0]),float(a[-1,0])],
                             'step_min_max_deg':[float(np.diff(a[:,0]).min()),float(np.diff(a[:,0]).max())],
                             'intensity_range':[float(a[:,1].min()),float(a[:,1].max())]})
        return samples,{'path':str(Path(path).resolve()),'sha256':digest(Path(path)), 'sheet':sheet,'samples':metadata,
                        'intensity_status':'workbook numerical values, not assumed independent Poisson counts'}
    finally:wb.close()


def profile(x,p,reference):
    height,center,fwhm,eta,b0,b1=p
    z=2*(x-center)/fwhm
    return b0+b1*(x-reference)+height*(eta/(1+z*z)+(1-eta)*np.exp(-np.log(2)*z*z))


def fit_peak(x,y,target,halfwidth=.35):
    m=abs(x-target)<=halfwidth;xx,yy=x[m],y[m]
    if len(xx)<12:raise ValueError('Insufficient points for a six-parameter local fit')
    bg=np.median(np.r_[yy[:3],yy[-3:]]);amplitude=max(yy.max()-bg,1)
    p0=[amplitude,float(xx[yy.argmax()]),.12,.5,float(bg),0]
    fit=least_squares(lambda p:profile(xx,p,target)-yy,p0,
                      bounds=([0,target-halfwidth/2,.025,0,-np.inf,-np.inf],
                              [np.inf,target+halfwidth/2,1.5,1,np.inf,np.inf]),
                      x_scale='jac',max_nfev=3000)
    if not fit.success:raise RuntimeError(f'Local profile fit failed: {fit.message}')
    residual=fit.fun;dof=len(xx)-6
    covariance=np.linalg.pinv(fit.jac.T@fit.jac)*np.dot(residual,residual)/dof
    return {'target_deg':target,'center_two_theta_deg':float(fit.x[1]),
            'center_formal_se_deg':float(np.sqrt(max(covariance[1,1],0))),
            'observed_envelope_fwhm_deg':float(fit.x[2]),'pseudo_voigt_eta':float(fit.x[3]),
            'height_workbook_units':float(fit.x[0]),'background_at_target':float(fit.x[4]),
            'background_slope':float(fit.x[5]),'halfwidth_deg':halfwidth,'points':len(xx),
            'residual_rms':float(np.sqrt(np.mean(residual**2))),
            'normalized_residual_rms':float(np.sqrt(np.mean(residual**2))/fit.x[0]),
            'jacobian_rank':int(np.linalg.matrix_rank(fit.jac)),
            'parameters':fit.x.tolist(),
            'uncertainty_note':'Formal local least-squares SE only; excludes calibration, processing correlation, peak overlap, Kalpha shape and model bias.'}


def sine_coordinate(two_theta):
    return 2*np.sin(np.deg2rad(np.asarray(two_theta)/2))


def delta_from_001_020_satellite(angles):
    # For unique b: (lambda/d_02+delta,1)^2 = (2+delta)^2 B + C.
    s001,s020,ssat=sine_coordinate(angles)
    ratio=(ssat**2-s001**2)/s020**2
    if ratio<=0:raise ValueError('Angles cannot support this conditional reflection assignment')
    return float(2*np.sqrt(ratio)-2)


def propagated_delta_se(angles,standard_errors):
    derivative=[]
    for i in range(3):
        d=np.zeros(3);d[i]=1e-5
        derivative.append((delta_from_001_020_satellite(np.array(angles)+d)-delta_from_001_020_satellite(np.array(angles)-d))/2e-5)
    return float(np.linalg.norm(np.array(derivative)*standard_errors))


def conditional_metric(centers):
    # Overdetermined A/LFSO metric in units of lambda. No physical lattice scale is claimed.
    labels=['anchor_001_assumed','anchor_020_assumed','feature_1','feature_2','feature_3','feature_4']
    hkls=np.array([[0,0,1],[0,2,0],[1,1,0],[-1,1,1],[0,2,1],[1,1,1]],float)
    design=np.c_[hkls[:,0]**2,hkls[:,1]**2,hkls[:,2]**2,hkls[:,0]*hkls[:,2]]
    observed=np.array([centers[k]['center_two_theta_deg'] for k in labels])
    metric=np.linalg.lstsq(design,sine_coordinate(observed)**2,rcond=None)[0]
    predicted=2*np.rad2deg(np.arcsin(np.sqrt(design@metric)/2))
    return {'labels':labels,'conditional_hkl':hkls.tolist(), 'metric_coefficients_scaled_by_lambda_squared':metric.tolist(),
            'metric_order':['h_squared','k_squared','l_squared','h_times_l'],
            'observed_two_theta_deg':observed.tolist(),'predicted_two_theta_deg':predicted.tolist(),
            'residual_two_theta_deg':(predicted-observed).tolist(),'rms_two_theta_deg':float(np.sqrt(np.mean((predicted-observed)**2))),
            'degrees_of_freedom':2,
            'scope':'Conditional indexing check only; no zero-shift correction or phase-purity assumption verified.'}


def run(config_path):
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve()
    out.mkdir(parents=True,exist_ok=False);write_json(out/'config.json',cfg)
    samples,manifest=read_workbook(cfg['workbook'],cfg['sheet']);write_json(out/'input_manifest.json',manifest)
    fits={};detections={};colors={'LFSO':'gray','LFNSO21-A':'crimson','LFNSO21-B':'navy'}
    fig,axes=plt.subplots(2,1,figsize=(12,8))
    for name,a in samples.items():
        x,y=a.T
        with (out/f'{name}_profile.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['two_theta_deg','intensity_as_supplied']);w.writerows(a)
        noise=1.4826*np.median(abs(np.diff(y)-np.median(np.diff(y))))/np.sqrt(2)
        prominence=max(6*noise,.005*(y.max()-y.min()))
        sm=savgol_filter(y,7,2);idx,properties=find_peaks(sm,prominence=prominence,distance=8)
        detections[name]={'noise_estimate':float(noise),'prominence_threshold':float(prominence),
                          'purpose':'discovery candidates only; fit uses unsmoothed data',
                          'two_theta_candidates_deg':x[idx].tolist()}
        fits[name]={}
        for label,target in cfg['peak_targets_deg'][name].items():
            variants=[fit_peak(x,y,target,h) for h in cfg['profile_halfwidths_deg']]
            fit=variants[1].copy();fit['window_center_min_max_deg']=[min(v['center_two_theta_deg'] for v in variants),max(v['center_two_theta_deg'] for v in variants)]
            fit['window_variants']=variants
            flags=[]
            if np.ptp(fit['window_center_min_max_deg'])>.02:flags.append('center changes by >0.02 deg with fit window')
            if fit['normalized_residual_rms']>.08:flags.append('profile residual exceeds 8 percent of fitted height')
            if fit['jacobian_rank']<6:flags.append('rank-deficient local fit')
            fit['quality_flags']=flags
            fit['position_use']='caution: overlap/shape mismatch possible' if flags else 'local envelope position; systematic uncertainty remains'
            wavelength=cfg.get('wavelength_A')
            if wavelength is not None:
                if wavelength<=0:raise ValueError('Positive wavelength required')
                fit['d_A']=float(wavelength/sine_coordinate(fit['center_two_theta_deg']))
                fit['g_nm_inv']=10/fit['d_A']
                fit['wavelength_assumption_status']=cfg['wavelength_source']
                fit['d_A_wavelength_sensitivity']={str(w):float(w/sine_coordinate(fit['center_two_theta_deg'])) for w in cfg.get('wavelength_sensitivity_A',[])}
            fits[name][label]=fit
        for ax in axes:ax.plot(x,y,label=name,color=colors[name],lw=.8)
    axes[0].set(xlim=(10,80),yscale='log',ylabel='Workbook intensity');axes[0].legend()
    axes[1].set(xlim=(17,36),ylim=(100,2800),xlabel='2theta (degrees)',ylabel='Workbook intensity')
    fig.tight_layout();fig.savefig(out/'profiles.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for ax,(name,label) in zip(axes.flat,[(n,l) for n in ['LFNSO21-A','LFNSO21-B'] for l in ['feature_1','feature_2','feature_3']]):
        fit=fits[name][label];x,y=samples[name].T;m=abs(x-fit['target_deg'])<=fit['halfwidth_deg']
        ax.plot(x[m],y[m],'.',ms=3,color=colors[name],label='data')
        ax.plot(x[m],profile(x[m],fit['parameters'],fit['target_deg']),'k-',lw=1,label='envelope fit')
        ax.set(title=f"{name}: {fit['center_two_theta_deg']:.4f} deg",xlabel='2theta (degrees)');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'superstructure_peak_fits.png',dpi=160);plt.close(fig)
    b=fits['LFNSO21-B'];keys=['anchor_001_assumed','anchor_020_assumed','feature_3']
    angles=[b[k]['center_two_theta_deg'] for k in keys];ses=[b[k]['center_formal_se_deg'] for k in keys]
    delta=delta_from_001_020_satellite(angles)
    window_deltas=[delta_from_001_020_satellite([b[k]['window_variants'][i]['center_two_theta_deg'] for k in keys]) for i in range(3)]
    s001,s020=sine_coordinate(angles[:2]);conditional={}
    for d in [0,1/6,.25,1/3,delta]:
        pred=2*np.rad2deg(np.arcsin(np.sqrt(s001**2+(2+d)**2*s020**2/4)/2))
        conditional[str(d)]={'predicted_feature3_two_theta_deg':float(pred),'residual_to_observed_deg':float(pred-angles[2])}
    shift={k:float(b[k]['center_two_theta_deg']-fits['LFNSO21-A'][k]['center_two_theta_deg']) for k in ['anchor_001_assumed','anchor_020_assumed','high_35','high_37']}
    metric_checks={name:conditional_metric(fits[name]) for name in ['LFSO','LFNSO21-A']}
    hypotheses=compare(b,metric_checks['LFNSO21-A']['metric_coefficients_scaled_by_lambda_squared'],cfg.get('wavelength_A'))
    zero_diagnostic=[{'subtracted_zero_shift_deg':z,'delta':delta_from_001_020_satellite(np.array(angles)-z)} for z in [-.05,0,.05]]
    absolute=None
    if cfg.get('wavelength_A'):
        bstar=5*sine_coordinate(angles[1])/cfg['wavelength_A']
        absolute={'bstar_from_assumed_020_nm_inv':float(bstar),'q_nm_inv':float(delta*bstar),'real_space_period_nm':float(1/(delta*bstar)),
                  'wavelength_A':cfg['wavelength_A'],'status':'Cu Kalpha assumption plus conditional hkl'}
    comparison=None
    ed_path=cfg.get('prior_position_result')
    if ed_path:
        prior_path=(cp.parent/ed_path).resolve();prior=json.loads(prior_path.read_text());ed=prior['ed_empirical_modulation']
        comparison={'source':str(prior_path),'sha256':digest(prior_path),'empirical_ED_delta':ed['median_delta_v'],
                    'empirical_ED_q_nm_inv':ed['q_norm_nm_inv'],'empirical_ED_v_nm_inv':ed['basis_norms_nm_inv'][1],
                    'source_basis_status':'v tentatively b*; zone convention not verified',
                    'XRD_delta_minus_ED_delta':delta-ed['median_delta_v']}
        if absolute:
            comparison.update({'XRD_q_minus_ED_q_nm_inv':absolute['q_nm_inv']-ed['q_norm_nm_inv'],
                               'XRD_bstar_over_ED_v':absolute['bstar_from_assumed_020_nm_inv']/ed['basis_norms_nm_inv'][1]})
    overlaps={}
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,name in zip(axes,['LFNSO21-A','LFNSO21-B']):
        x,y=samples[name].T;overlaps[name]=fit_group(x,y)
        grp=overlaps[name];lo,hi=grp['region_deg'];mask=(x>=lo)&(x<=hi)
        ax.plot(x[mask],y[mask],'.',color=colors[name],label='data')
        for model in grp['models']:
            ax.plot(x[mask],group_profile(x[mask],model['parameters'],model['component_count'],grp['reference_deg']),label=str(model['component_count'])+' envelope(s)')
        ax.set(title=name+' / 35-degree cluster',xlabel='2theta (degrees)');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'overlap_diagnostic.png',dpi=160);plt.close(fig)
    result={'overlapping_envelope_diagnostic':overlaps,'conditional_hypothesis_comparison':hypotheses,'absolute_q_if_assumed_wavelength':absolute,
            'ED_comparison':comparison,'zero_shift_sensitivity_diagnostic_not_uncertainty':zero_diagnostic,'input':manifest,'wavelength_A':cfg.get('wavelength_A'),'wavelength_source':cfg['wavelength_source'],
            'preprocessing_history':cfg['preprocessing_history'],'fits':fits,'automatic_peak_candidates':detections,
            'A_reference_indexing_checks':metric_checks,
            'B_minus_A_empirical_feature_shifts_deg':shift,
            'conditional_pure_bstar_modulation':{'delta':delta,'formal_propagated_se':propagated_delta_se(angles,ses),
                'window_variation_deltas':window_deltas,'formula':'delta = 2 sqrt((s_sat^2 - s_001^2)/s_020^2) - 2; s=2sin(theta)',
                'wavelength_required':False,'zero_shift_applied_deg':0,'zero_shift_status':'unverified',
                'assumed_reflections':['001','020','0,(2+delta),1'],
                'predictions':conditional,
                'limitation':'Conditional reflection assignment and pure b* modulation; not a unique structure or global 3D modulation refinement.'},
            'status':'raw-profile peak fitting and conditional metric checks; atomic structure not established',
            'gate1':'not_evaluated','gate2':'not_evaluated','stem_simulation':'not_run'}
    write_json(out/'result.json',result)
    with (out/'peak_centers.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['sample','feature','two_theta_deg','formal_se_deg','window_center_min_deg','window_center_max_deg','observed_fwhm_deg','d_A_at_provisional_wavelength'])
        for name,fs in fits.items():
            for label,r in fs.items():w.writerow([name,label,r['center_two_theta_deg'],r['center_formal_se_deg'],*r['window_center_min_max_deg'],r['observed_envelope_fwhm_deg'],r.get('d_A','')])
    body='<!doctype html><html lang="ko"><meta charset="utf-8"><title>LFNSO XRD numerical profiles</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:20px}pre{white-space:pre-wrap}img{max-width:100%}</style><h1>원본 엑셀 XRD 피크 위치 분석</h1><p>공유된 엑셀 수치에 단일 pseudo-Voigt envelope와 국소 선형 배경을 피팅했습니다. Profile fitting에는 평활화하지 않은 값을 사용했습니다. 정밀 기기 보정·Kα doublet 및 처리 이력은 별도로 확인해야 합니다. Formal SE는 이 모든 계통 오차를 포함하지 않습니다.</p><p>Cu Kα를 사용자 제안에 따라 가정했으며 1.5406 Å를 잠정 기준으로 씁니다. 단색화나 Kα₂ 제거 여부는 확인되지 않았습니다. 설정에서 파장을 지정하지 않으면 d-spacing을 출력하지 않습니다. 조건부 δ는 sin²θ 비율로 산출하여 파장에 의존하지 않습니다. 이는 반사 대응과 pure b* 가정에 의존하며 전체 구조 판정이 아닙니다.</p>'
    for p in sorted(out.glob('*.png')):body+='<img src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'">'
    summary={k:v for k,v in result.items() if k not in ['fits','automatic_peak_candidates']}
    body+='<h2>가정 및 계산 결과</h2><pre>'+html.escape(json.dumps(summary,indent=2,ensure_ascii=False))+'</pre><h2>피크 위치</h2><pre>'+html.escape((out/'peak_centers.csv').read_text())+'</pre></html>'
    (out/'report.html').write_text(body)
    print(json.dumps({'output':str(out),'samples':len(samples),'fitted_features':sum(len(v) for v in fits.values()),'wavelength_A':cfg.get('wavelength_A'),'conditional_delta':delta}))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
