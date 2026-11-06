"""Conserved slab-site permutations with cached X-ray structure factors.

The only mixed site species is Fe0.67Ni0.33. No experimental intensity is
fabricated; simulated profiles are optimization predictions, not measured data.
"""
import argparse,base64,hashlib,html,json,time
from pathlib import Path
import numpy as np
from pymatgen.core import Structure,Element
from pymatgen.analysis.diffraction.xrd import XRDCalculator,ATOMIC_SCATTERING_PARAMS
from .model_generator import make_parent,validate,FE_NI
from .pipeline import write_json,digest


def form_factor(symbol,g):
    c=np.asarray(ATOMIC_SCATTERING_PARAMS[symbol]);s2=(np.asarray(g)/2)**2
    return Element(symbol).Z-41.78214*s2*np.sum(c[:,0]*np.exp(-c[:,1]*s2[:,None]),axis=1)


class DiffractionCache:
    """Same f(s), LP and Friedel factors as pymatgen; updated only at swapped sites."""
    def __init__(self,structure,wavelength=1.5406,angle_range=(16,46)):
        self.structure=structure.copy();self.wavelength=wavelength
        reciprocal=structure.lattice.reciprocal_lattice_crystallographic
        lower,upper=2*np.sin(np.deg2rad(np.array(angle_range)/2))/wavelength
        pts=reciprocal.get_points_in_sphere([[0,0,0]],[0,0,0],upper)
        pairs=[(np.rint(p[0]).astype(int),p[1]) for p in pts if p[1]>=lower and tuple(np.rint(p[0]).astype(int))>(0,0,0)]
        pairs.sort(key=lambda p:(p[1],*-p[0]))
        self.hkl=np.array([p[0] for p in pairs]);g=np.array([p[1] for p in pairs])
        theta=np.arcsin(wavelength*g/2);angles=np.rad2deg(theta)*2
        keys=np.rint(angles/XRDCalculator.TWO_THETA_TOL).astype(np.int64)
        _,first,self.groups=np.unique(keys,return_index=True,return_inverse=True)
        self.angles=angles[first]
        self.lp=2*(1+np.cos(2*theta)**2)/(np.sin(theta)**2*np.cos(theta))
        self.slab=np.array([i for i,t in enumerate(structure) if abs(t.frac_coords[2]-np.rint(t.frac_coords[2]))<1e-7 and t.species_string!='O'])
        if len(self.slab)%3:raise ValueError('Slab sites must support equal Li/Sb/Fe-host populations')
        phase=np.exp(2j*np.pi*(self.hkl@structure.frac_coords.T))
        self.phase=phase[:,self.slab]
        fs={e:form_factor(e,g) for e in ['Li','O','Sb','Fe','Ni']}
        self.factors=np.array([fs['Li'],fs['Sb'],.67*fs['Fe']+.33*fs['Ni']])
        self.constant=np.zeros(len(g),complex)
        fixed=set(range(len(structure)))-set(self.slab)
        for i in fixed:
            for element,amount in structure[i].species.items():self.constant+=phase[:,i]*fs[element.symbol]*amount
        self.initial=np.array([0 if structure[i].species_string=='Li' else 1 if structure[i].species_string=='Sb' else 2 for i in self.slab])

    def factors_for(self,labels):
        return self.constant+np.sum(self.phase*self.factors[np.asarray(labels)].T,axis=1)

    def swap(self,f,labels,i,j):
        return f+(self.phase[:,i]-self.phase[:,j])*(self.factors[labels[j]]-self.factors[labels[i]])

    def intensities(self,f):
        return np.bincount(self.groups,weights=np.abs(f)**2*self.lp,minlength=len(self.angles))

    def to_structure(self,labels):
        s=self.structure.copy()
        for i,label in zip(self.slab,labels):
            s.replace(int(i),['Li','Sb',FE_NI][int(label)],properties=s[i].properties.copy())
        for t in s:t.label=t.species_string
        return s

    def verify_against_pymatgen(self,labels):
        p=XRDCalculator(wavelength=self.wavelength).get_pattern(self.to_structure(labels),two_theta_range=(16,46))
        intensity=self.intensities(self.factors_for(labels));intensity=intensity/intensity.max()*100
        mask=intensity>XRDCalculator.SCALED_INTENSITY_TOL
        if len(p.x)!=sum(mask) or not np.allclose(self.angles[mask],p.x,atol=1e-7) or not np.allclose(intensity[mask],p.y,rtol=1e-6,atol=1e-6):
            raise AssertionError('Cached diffraction failed independent pymatgen comparison')
        return {'peak_count':len(p.x),'max_relative_intensity_difference':float(np.max(abs(intensity[mask]-p.y)))}


class ProfileObjective:
    def __init__(self,angles,x,y,width):
        self.x=np.asarray(x);self.y=np.asarray(y);self.width=width
        self.background=np.c_[np.ones(len(x)),x-np.mean(x)]
        q,_=np.linalg.qr(self.background)
        self.kernel=np.exp(-4*np.log(2)*((x[:,None]-angles[None,:])/width)**2)/width
        residual_kernel=self.kernel-q@(q.T@self.kernel)
        self.yr=y-q@(q.T@y);self.sst=float(self.yr@self.yr)
        self.linear=residual_kernel.T@self.yr
        self.quad=residual_kernel.T@residual_kernel

    def correlation(self,intensity):
        variance=max(float(intensity@self.quad@intensity),1e-30)
        return float(np.clip((self.linear@intensity)/np.sqrt(variance*self.sst),-1,1))

    def evaluate(self,intensity):
        model=self.kernel@intensity
        projected=model-self.background@np.linalg.lstsq(self.background,model,rcond=None)[0]
        scale=max(0,float(projected@self.yr)/max(float(projected@projected),1e-30))
        bg=np.linalg.lstsq(self.background,self.y-scale*model,rcond=None)[0]
        fit=scale*model+self.background@bg;residual=fit-self.y
        return {'signed_correlation':self.correlation(intensity),'positive_scale':scale,'linear_background':bg.tolist(),
                'R2_vs_linear_background':float(1-np.dot(residual,residual)/self.sst),
                'normalized_RMS_by_data_range':float(np.sqrt(np.mean(residual**2))/np.ptp(self.y)),
                'profile':fit.tolist(),'width_deg':self.width}


def chain(cache,objective,seed,steps,path,signature):
    rng=np.random.default_rng(seed);n=len(cache.slab);history=[]
    if path.exists():
        checkpoint=json.loads(path.read_text())
        if checkpoint['signature']!=signature:raise ValueError('Checkpoint does not match data/config/code')
        labels=np.array(checkpoint['current_labels']);best=np.array(checkpoint['best_labels']);start=checkpoint['step'];rng.bit_generator.state=checkpoint['rng_state'];history=checkpoint['history']
    else:
        labels=rng.permutation(np.repeat(np.arange(3),n//3));best=labels.copy();start=0
    f=cache.factors_for(labels);score=objective.correlation(cache.intensities(f));bestscore=objective.correlation(cache.intensities(cache.factors_for(best)))
    accepted=0;new_scores=0
    for step in range(start,steps):
        i,j=rng.choice(n,2,replace=False)
        if labels[i]!=labels[j]:
            trial=cache.swap(f,labels,i,j);value=objective.correlation(cache.intensities(trial));new_scores+=1
            temperature=.035*(.0005/.035)**(step/max(steps-1,1))
            if value>=score or rng.random()<np.exp((value-score)/temperature):
                labels[i],labels[j]=labels[j],labels[i];f=trial;score=value;accepted+=1
                if score>bestscore+1e-12:bestscore=score;best=labels.copy();history.append({'step':step+1,'best_correlation':bestscore})
        if (step+1)%250==0 or step==steps-1:
            # Rebuild to prevent accumulated incremental round-off, and persist restartable state.
            f=cache.factors_for(labels);score=objective.correlation(cache.intensities(f))
            write_json(path,{'signature':signature,'step':step+1,'current_labels':labels.tolist(),'best_labels':best.tolist(),
                             'rng_state':rng.bit_generator.state,'history':history,'seed':seed,'best_correlation':bestscore})
    return best,{'proposal_steps':steps,'new_diffraction_evaluations_this_invocation':new_scores,'accepted_swaps_this_invocation':accepted,'history':history}


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=True)
    source=Path(__file__)
    fingerprint={'config':cfg,'code_sha256':digest(source),'raw_profile_sha256':digest(Path(cfg['experimental_csv'])),'template_sha256':digest(Path(cfg['template']))}
    signature=hashlib.sha256(json.dumps(fingerprint,sort_keys=True).encode()).hexdigest()
    if (out/'provenance.json').exists() and json.loads((out/'provenance.json').read_text())['signature']!=signature:raise ValueError('Use a new output directory for changed code/config/data')
    write_json(out/'provenance.json',{'signature':signature,**fingerprint})
    data=np.loadtxt(cfg['experimental_csv'],delimiter=',',skiprows=1)
    template=Structure.from_file(cfg['template']);rows=[];t0=time.monotonic()
    for metric_name,metric in cfg['metrics'].items():
        for period in cfg['periods']:
            parent=make_parent(template,[8,9],metric);parent.make_supercell([1,period,1]);cache=DiffractionCache(parent,cfg['wavelength_A'])
            masks={name:(data[:,0]>=bounds[0])&(data[:,0]<=bounds[1]) for name,bounds in cfg['windows'].items()}
            objectives={name:ProfileObjective(cache.angles,*data[mask].T,cfg['width_deg']) for name,mask in masks.items()}
            print('cache',metric_name,period,'sites',len(parent),'reflections',len(cache.hkl),'peak_groups',len(cache.angles),flush=True)
            verification=cache.verify_against_pymatgen(cache.initial)
            for seed in cfg['seeds']:
                name=f'{metric_name}_{period}b_seed{seed}';folder=out/name;folder.mkdir(exist_ok=True)
                labels,log=chain(cache,objectives['train'],seed,cfg['steps'],folder/'checkpoint.json',signature)
                f=cache.factors_for(labels);intensities=cache.intensities(f);s=cache.to_structure(labels)
                independent=cache.verify_against_pymatgen(labels)
                checks=validate(s,2*period)
                s.to(filename=str(folder/'model.cif'),fmt='cif',symprec=None)
                scores={name:obj.evaluate(intensities) for name,obj in objectives.items()}
                baseline={name:obj.evaluate(cache.intensities(cache.factors_for(cache.initial))) for name,obj in objectives.items()}
                # Each held-out fit is evaluated once after training; its scale/background are separate nuisance fits.
                serial_scores={name:{k:v for k,v in score.items() if k!='profile'} for name,score in scores.items()}
                write_json(folder/'xrd.json',{'two_theta_deg':cache.angles.tolist(),'intensity_arbitrary_common_units':intensities.tolist(),'wavelength_A':cfg['wavelength_A']})
                row={'id':name,'metric_name':metric_name,'metric_parameters':metric,'period_b':period,'seed':seed,
                     'labels_Li_Sb_FeNi':labels.tolist(),'slab_indices':cache.slab.tolist(),'scores':serial_scores,
                     'initial_ordered_control_scores':{name:{k:v for k,v in score.items() if k!='profile'} for name,score in baseline.items()},
                     'validation':checks,'independent_diffraction_validation':independent,'search_log':log,
                     'hypothesis':'Conserved permutations across six octahedral slab sites per parent cell; interlayer Li and O fixed, Fe/Ni share each Fe-host site.',
                     'sampled_configuration_status':'not an exhaustive or proven global optimum; no symmetry-distinct count claimed',
                     'gate1':'not_evaluated','gate2':'not_evaluated','stem':'not_run'}
                write_json(folder/'metadata.json',row);rows.append(row)
                fig,axes=plt.subplots(1,len(scores),figsize=(15,4))
                for ax,(window,score) in zip(axes,scores.items()):
                    obj=objectives[window];ax.plot(obj.x,obj.y,'k-',lw=.7,label='B data');ax.plot(obj.x,score['profile'],'r-',lw=.8,label='predicted profile + baseline')
                    ax.set(title=f"{window}: R²={score['R2_vs_linear_background']:.3f}",xlabel='2theta (degrees)');ax.legend(fontsize=7)
                fig.suptitle(name+' (holdouts not optimized)');fig.tight_layout();fig.savefig(folder/'comparison.png',dpi=130);plt.close(fig)
                print(name,{k:round(v['R2_vs_linear_background'],4) for k,v in scores.items()},flush=True)
    summary={'config':cfg,'elapsed_seconds':time.monotonic()-t0,'provenance_signature':signature,
             'rows':[{k:v for k,v in row.items() if k not in ['validation','search_log','labels_Li_Sb_FeNi','slab_indices']}|{'geometry_summary':{k:v for k,v in row['validation'].items() if k!='coordination'}} for row in rows],
             'limitations':['Raw composition is exact in average Fe/Ni occupancy, not integer Fe/Ni local arrangements.',
                            'Metric options and train-window choices use earlier exploratory XRD analysis; holdouts were not used in this optimization.',
                            'Gaussian width is fixed, wavelength provisional, zero shift unrefined, no energy relaxation.',
                            'Local configuration space is large; multiple seeds do not establish completeness.',
                            'Per-window positive scale and linear baseline are fitted separately; holdout score tests shape, not transfer of absolute scale.',
                            'No dynamical ED or STEM validation; model probabilities and structural uniqueness are not established.']}
    write_json(out/'summary.json',summary)
    page='<!doctype html><meta charset="utf-8"><title>LFNSO constrained occupancy search</title><style>body{font:16px system-ui;max-width:1300px;margin:25px auto}img{max-width:100%}pre{white-space:pre-wrap}</style><h1>Slab-site permutation search</h1><p>Fixed O/interlayer coordinates, exact average composition, Fe-site Ni co-occupancy. These are trial structures, not established models.</p>'
    for row in sorted(rows,key=lambda r:r['scores']['train']['R2_vs_linear_background'],reverse=True):
        page+='<h2>'+row['id']+'</h2><img src="data:image/png;base64,'+base64.b64encode((out/row['id']/'comparison.png').read_bytes()).decode()+'">'
    page+='<pre>'+html.escape(json.dumps(summary,indent=2,ensure_ascii=False))+'</pre>';(out/'report.html').write_text(page)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
