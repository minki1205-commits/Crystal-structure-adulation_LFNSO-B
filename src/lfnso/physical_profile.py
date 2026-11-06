"""Positive structural scale and nonnegative linear-background endpoint fitting."""
import numpy as np
from scipy.optimize import nnls
from .occupational_search import ProfileObjective


class NonnegativeBaseline:
    def __init__(self,x,y):
        self.x=np.asarray(x);self.y=np.asarray(y)
        t=(self.x-self.x.min())/(self.x.max()-self.x.min())
        self.background=np.c_[1-t,t]
        bg,_=nnls(self.background,self.y);residual=self.y-self.background@bg
        self.sst=float(residual@residual)
    def fit(self,profile,details=False):
        profile=np.asarray(profile);normalizer=max(float(np.linalg.norm(profile)),1e-30)
        design=np.c_[profile/normalizer,self.background]
        coef,rnorm=nnls(design,self.y)
        r2=float(1-rnorm*rnorm/self.sst)
        if not details:return max(0.,r2)**.5
        fitted=design@coef
        return {'R2_vs_nonnegative_linear_background':r2,'positive_scale':float(coef[0]/normalizer),
                'background_endpoint_values':coef[1:].tolist(),'background_nonnegative':True,
                'normalized_RMS_by_data_range':float(rnorm/np.sqrt(len(profile))/np.ptp(self.y)),
                'profile':fitted.tolist()}


class NonnegativeProfileObjective:
    def __init__(self,angles,x,y,width):
        self.unconstrained=ProfileObjective(angles,x,y,width)
        self.baseline=NonnegativeBaseline(x,y);self.width=width
    def evaluate(self,intensity):
        result=self.baseline.fit(self.unconstrained.kernel@intensity,details=True)
        result['width_deg']=self.width
        return result
