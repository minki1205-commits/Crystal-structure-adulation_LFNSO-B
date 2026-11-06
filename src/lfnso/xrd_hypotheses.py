"""Conditional reciprocal-metric diagnostics, not chemical-structure refinement."""
import numpy as np
from scipy.optimize import least_squares


def design(hkl):
    h=np.asarray(hkl,float)
    return np.c_[h[:,0]**2,h[:,1]**2,h[:,2]**2,h[:,0]*h[:,2]]


def angles_from_metric(hkl,metric):
    squared=design(hkl)@np.asarray(metric)
    if (squared<=0).any() or (squared>=4).any():raise ValueError('Reflections outside valid Bragg range')
    return 2*np.rad2deg(np.arcsin(np.sqrt(squared)/2))


def physical_metric(metric,wavelength):
    aa,bb,cc,dd=metric
    if min(aa,bb,cc)<=0 or dd*dd>=4*aa*cc:raise ValueError('Metric is not positive definite')
    beta=np.arccos(-dd/(2*np.sqrt(aa*cc)))
    return {'a_A':float(wavelength/(np.sqrt(aa)*np.sin(beta))),'b_A':float(wavelength/np.sqrt(bb)),
            'c_A':float(wavelength/(np.sqrt(cc)*np.sin(beta))),'beta_deg':float(np.rad2deg(beta)),
            'status':'algebraic conditional solution, not refined physical structure'}


def compare(fits,a_metric,wavelength=None):
    labels=['anchor_001_assumed','anchor_020_assumed','feature_1','feature_2','feature_3']
    hkl=np.array([[0,0,1],[0,2,0],[1,1,0],[-1,1,1],[0,2,1]],float)
    observed=np.array([fits[k]['center_two_theta_deg'] for k in labels])
    sine=2*np.sin(np.deg2rad(observed/2));squared=sine**2
    h0=np.linalg.lstsq(design(hkl),squared,rcond=None)[0]
    p0=angles_from_metric(hkl,h0)
    def pred(delta,scale=1):
        h=hkl.copy();h[2:,1]+=delta
        return angles_from_metric(h,np.asarray(a_metric)/scale**2)
    fixed=[]
    for delta in [1/6,.25,1/3]:
        p=pred(delta)
        fixed.append({'delta':delta,'predicted_angles_deg':p.tolist(),'residual_deg':(p-observed).tolist(),
                      'rms_all_five_deg':float(np.sqrt(np.mean((p-observed)**2))),
                      'rms_three_features_deg':float(np.sqrt(np.mean((p[2:]-observed[2:])**2)))})
    opt=least_squares(lambda z:pred(z[0],z[1])-observed,[.28,1],bounds=([0,.97],[.5,1.03]))
    delta=2*np.sqrt((squared[4]-squared[0])/squared[1])-2
    bb=squared[1]/4;cc=squared[0];aa=squared[2]-bb*(1+delta)**2;dd=squared[2]+cc-squared[3]
    coupled=[aa,bb,cc,dd]
    algebraic={'delta':float(delta),'metric_lambda_squared':list(map(float,coupled)),
               'fit_targets':labels,'observations':5,'free_parameters':5,'residual_degrees_of_freedom':0,
               'warning':'Exact interpolation of five assigned peak positions is not validation.'}
    if wavelength is not None:algebraic['cell_if_wavelength_assumed']=physical_metric(coupled,wavelength)
    # Explicit held-out competing reflection families, without assigning a blended peak to a single hkl.
    heldout_hkl=[[1,3,0],[-2,0,1],[2,0,0],[0,0,2]]
    algebraic['held_out_predictions']=[{'parent_hkl':h,'two_theta_deg':float(t)} for h,t in zip(heldout_hkl,angles_from_metric(heldout_hkl,coupled))]
    return {'observed_labels':labels,'observed_angles_deg':observed.tolist(),
            'H0_same_hkl_metric_only':{'metric_lambda_squared':h0.tolist(),'predicted_angles_deg':p0.tolist(),
                'residual_deg':(p0-observed).tolist(),'rms_deg':float(np.sqrt(np.mean((p0-observed)**2))),
                'fit_space':'unweighted linear fit in 4 sin²(theta)','scope':'Only the specified same-hkl mapping is tested.'},
            'fixed_A_metric_satellites':fixed,
            'A_metric_free_delta_and_uniform_scale':{'delta':float(opt.x[0]),'length_scale':float(opt.x[1]),
                'residual_deg':opt.fun.tolist(),'rms_deg':float(np.sqrt(np.mean(opt.fun**2))),
                'free_parameters':2,'observations':5,'scale_bounds':[.97,1.03]},
            'B_free_metric_plus_delta':algebraic,
            'atomic_structure_ranking':'not performed; missing chemistry/site assignments and intensity model'}
