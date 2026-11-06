"""Empirical decomposition of overlapping envelopes, not phase/Kalpha assignment."""
import numpy as np
from scipy.optimize import least_squares


def group_profile(x,p,components,reference):
    value=p[-2]+p[-1]*(x-reference)
    for i in range(components):
        height,center,width,eta=p[4*i:4*i+4]
        z=2*(x-center)/width
        value=value+height*(eta/(1+z*z)+(1-eta)*np.exp(-np.log(2)*z*z))
    return value


def fit_group(x,y,region=(34.85,35.75)):
    m=(x>=region[0])&(x<=region[1]);xx,yy=x[m],y[m];reference=float(np.mean(region))
    if len(xx)<25:raise ValueError('Too few points for overlap comparison')
    bg=float(np.median(np.r_[yy[:3],yy[-3:]]));height=max(float(yy.max()-bg),1)
    results=[]
    for count in [1,2]:
        guesses=[[float(xx[yy.argmax()])]] if count==1 else [[35.12,35.36],[35.20,35.43],[35.28,35.40]]
        solutions=[]
        for centers in guesses:
            lower=[];upper=[];initial=[]
            for c in centers:
                initial.extend([height/count,c,.13,.5])
                lower.extend([0,region[0],.025,0]);upper.extend([np.inf,region[1],.7,1])
            initial.extend([bg,0]);lower.extend([-np.inf,-np.inf]);upper.extend([np.inf,np.inf])
            opt=least_squares(lambda p:group_profile(xx,p,count,reference)-yy,initial,bounds=(lower,upper),x_scale='jac',max_nfev=5000)
            if opt.success:solutions.append(opt)
        if not solutions:raise RuntimeError('Overlap envelope fits did not converge')
        best=min(solutions,key=lambda p:np.dot(p.fun,p.fun));rss=float(np.dot(best.fun,best.fun));n=len(xx);k=4*count+2
        aicc=n*np.log(max(rss/n,np.finfo(float).tiny))+2*k+2*k*(k+1)/(n-k-1)
        fitted=[{'center_deg':float(best.x[4*i+1]),'height':float(best.x[4*i]),'fwhm_deg':float(best.x[4*i+2]),'eta':float(best.x[4*i+3])} for i in range(count)]
        results.append({'component_count':count,'components':sorted(fitted,key=lambda r:r['center_deg']),
                        'rss':rss,'aicc_conditional_gaussian_residual':float(aicc),'parameters':best.x.tolist(),
                        'residual_rms':float(np.sqrt(rss/n)),'parameter_count':k})
    return {'region_deg':list(region),'points':len(xx),'reference_deg':reference,'models':results,
            'aicc_two_minus_one':results[1]['aicc_conditional_gaussian_residual']-results[0]['aicc_conditional_gaussian_residual'],
            'meaning':'Empirical envelope decomposition only. A lower AICc does not identify phases, hkl, Kalpha doublets or a unique component count.',
            'selection_note':'Region and multiple start centers chosen for the observed 35-degree cluster; not blind whole-pattern analysis.'}
