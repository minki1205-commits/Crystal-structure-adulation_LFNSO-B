"""Fixed-coordinate 3D occupational exploration with explicit XRD/ED tradeoff.

ED objective uses bounded geometric support, not measured electron intensities.
Within-layer swaps preserve each layer's Li:Sb:Fe-host counts.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from pymatgen.core import Lattice
from pymatgen.io.cif import CifParser
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from .occupational_search import form_factor,ProfileObjective
from .coupled_search import powder,verify_powder
from .coupled_models import plane_reflections_3d
from .ed_matching import electron_factor,transform,in_view,match_one_to_one,electron_intensity
from .ed_metric_diagnostic import apply_correction
from .model_generator import FE_NI,validate
from .pipeline import write_json,digest
from .physical_profile import NonnegativeBaseline,NonnegativeProfileObjective


class SwapFactors:
    def __init__(self,structure,hk,g,slab,factor_function):
        self.structure=structure;self.slab=np.asarray(slab,int)
        phase=np.exp(2j*np.pi*(hk@structure.frac_coords.T));self.phase=phase[:,slab]
        fs={e:factor_function(e,g) for e in ['Li','Sb','Fe','Ni','O']}
        self.factors=np.array([fs['Li'],fs['Sb'],.67*fs['Fe']+.33*fs['Ni']])
        self.constant=np.zeros(len(hk),complex)
        for i in sorted(set(range(len(structure)))-set(slab)):
            for element,amount in structure[i].species.items():self.constant+=phase[:,i]*fs[element.symbol]*amount
    def rebuild(self,labels):return self.constant+np.sum(self.phase*self.factors[labels].T,axis=1)
    def swap(self,f,labels,i,j):return f+(self.phase[:,i]-self.phase[:,j])*(self.factors[labels[j]]-self.factors[labels[i]])


class JointCache:
    def __init__(self,structure,repeats,parent_lattice,dataset,registration,correction,x,y,width,wavelength,nonnegative_background=False):
        self.structure=structure;self.repeats=repeats;nz=repeats[2]
        self.nonnegative_background=nonnegative_background;self.nn_background=NonnegativeBaseline(x,y)
        self.slab=np.array([i for i,t in enumerate(structure) if t.species_string!='O' and abs(t.frac_coords[2]*nz-np.rint(t.frac_coords[2]*nz))<1e-7])
        self.labels=np.array([0 if structure[i].species_string=='Li' else 1 if 'Sb' in structure[i].species else 2 for i in self.slab])
        self.layers=np.array([int(round(structure[i].frac_coords[2]*nz))%nz for i in self.slab])
        self.layer_slots=[np.flatnonzero(self.layers==k) for k in range(nz)]
        lo,hi=2*np.sin(np.deg2rad(np.array([16,46])/2))/wavelength
        pts=structure.lattice.reciprocal_lattice_crystallographic.get_points_in_sphere([[0,0,0]],[0,0,0],hi)
        pairs=[(np.rint(p[0]).astype(int),p[1]) for p in pts if p[1]>=lo and tuple(np.rint(p[0]).astype(int))>(0,0,0)]
        pairs.sort(key=lambda p:(p[1],*-p[0]));hk=np.array([p[0] for p in pairs]);g=np.array([p[1] for p in pairs]);theta=np.arcsin(wavelength*g/2);angles=np.rad2deg(theta)*2
        _,first,self.groups=np.unique(np.rint(angles/XRDCalculator.TWO_THETA_TOL).astype(np.int64),return_index=True,return_inverse=True)
        self.angles=angles[first];self.lp=2*(1+np.cos(2*theta)**2)/(np.sin(theta)**2*np.cos(theta))
        self.xrd=SwapFactors(structure,hk,g,self.slab,form_factor)
        dense=np.exp(-4*np.log(2)*((x[:,None]-self.angles)/width)**2)/width
        dense[abs(x[:,None]-self.angles)>1.5]=0 # negligible < 5e-44 at 0.25 degree FWHM
        self.kernel=csr_matrix(dense);self.background=np.c_[np.ones(len(x)),x-x.mean()];q,_=np.linalg.qr(self.background);self.q=q;self.yr=y-q@(q.T@y);self.ynorm=np.linalg.norm(self.yr)
        eh,_,eg,xy=plane_reflections_3d(parent_lattice,repeats,[1,1,0],14)
        pred=apply_correction(transform(xy,registration['parameters'],registration['mirror']),correction);visible=in_view(pred,dataset)
        self.ed=SwapFactors(structure,eh[visible],np.linalg.norm(eg[visible],axis=1),self.slab,electron_factor)
        assignment=match_one_to_one(dataset['observed'],pred[visible],.08)
        self.support_ids=np.array([r['predicted_index'] for r in assignment['pairs']],int);self.nobs=len(dataset['observed']);self.npred=visible.sum()
        self.ed_positions=pred[visible];self.assignment=assignment
    def intensities(self,f):return np.bincount(self.groups,weights=abs(f)**2*self.lp,minlength=len(self.angles))
    def score(self,fx,fe,weight):
        px=self.kernel@self.intensities(fx)
        if self.nonnegative_background:correlation=self.nn_background.fit(px)
        else:
            px-=self.q@(self.q.T@px)
            correlation=float(px@self.yr/max(np.linalg.norm(px)*self.ynorm,1e-30))
        intensity=abs(fe)**2;soft=intensity/(intensity+.01*max(float(intensity.max()),1e-30))
        edsoft=float(2*soft[self.support_ids].sum()/(self.nobs+soft.sum()))
        return correlation+weight*edsoft,correlation,edsoft
    def to_structure(self,labels):
        s=self.structure.copy()
        for i,label in zip(self.slab,labels):s.replace(int(i),['Li','Sb',FE_NI][int(label)])
        for site in s:site.label=site.species_string
        return s


def search(cache,seed,steps,weight,checkpoint,signature):
    rng=np.random.default_rng(seed);labels=cache.labels.copy();best=labels.copy();history=[];start=0;evaluations=0
    if checkpoint.exists():
        c=json.loads(checkpoint.read_text())
        if c['signature']!=signature:raise ValueError('Checkpoint signature differs')
        labels=np.array(c['labels']);best=np.array(c['best']);rng.bit_generator.state=c['rng'];start=c['step'];history=c['history'];evaluations=c['evaluations']
    fx=cache.xrd.rebuild(labels);fe=cache.ed.rebuild(labels);score=cache.score(fx,fe,weight)[0];bestscore=cache.score(cache.xrd.rebuild(best),cache.ed.rebuild(best),weight)[0]
    for step in range(start,steps):
        slots=cache.layer_slots[rng.integers(len(cache.layer_slots))];i,j=rng.choice(slots,2,replace=False)
        if labels[i]!=labels[j]:
            tx=cache.xrd.swap(fx,labels,i,j);te=cache.ed.swap(fe,labels,i,j);v=cache.score(tx,te,weight);evaluations+=1
            temp=.015*(.0002/.015)**(step/max(steps-1,1))
            if v[0]>=score or rng.random()<np.exp((v[0]-score)/temp):
                labels[i],labels[j]=labels[j],labels[i];fx,fe,score=tx,te,v[0]
                if score>bestscore+1e-10:bestscore=score;best=labels.copy();history.append({'step':step+1,'objective':v[0],'XRD_correlation':v[1],'ED_soft_support':v[2]})
        if (step+1)%500==0 or step==steps-1:
            fx=cache.xrd.rebuild(labels);fe=cache.ed.rebuild(labels);score=cache.score(fx,fe,weight)[0]
            write_json(checkpoint,{'signature':signature,'step':step+1,'labels':labels.tolist(),'best':best.tolist(),'rng':rng.bit_generator.state,'history':history,'evaluations':evaluations})
    return best,{'proposal_steps':steps,'nontrivial_evaluations':evaluations,'history':history}


def run(config_path):
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=True)
    edp=(cp.parent/cfg['ed_summary']).resolve();afp=(cp.parent/cfg['affine_summary']).resolve();ed=json.loads(edp.read_text());affine=json.loads(afp.read_text());dataset=next(d for d in ed['datasets'] if d['id']=='dp3_110')
    raw=np.loadtxt(cfg['xrd_csv'],delimiter=',',skiprows=1);mask=(raw[:,0]>=20.4)&(raw[:,0]<=34.6)
    provenance={'config':cfg,'source_sha256':{p.name:digest(p) for p in [Path(__file__),Path(__file__).with_name('coupled_search.py'),Path(__file__).with_name('ed_matching.py'),Path(__file__).with_name('coupled_models.py'),Path(__file__).with_name('physical_profile.py')]},'ed_sha256':digest(edp),'affine_sha256':digest(afp),'xrd_sha256':digest(Path(cfg['xrd_csv'])),'base_hashes':{b['id']:digest((cp.parent/b['cif']).resolve()) for b in cfg['bases']}}
    signature=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
    if (out/'provenance.json').exists() and json.loads((out/'provenance.json').read_text())['signature']!=signature:raise ValueError('Use a new output directory')
    write_json(out/'provenance.json',{'signature':signature,**provenance});rows=[];t0=time.monotonic()
    for base_entry in cfg['bases']:
        base=CifParser((cp.parent/base_entry['cif']).resolve(),frac_tolerance=0).parse_structures(primitive=False)[0];base.make_supercell([1,1,cfg['c_repeat']]);repeats=[1,base_entry['b_repeat'],cfg['c_repeat']]
        lattice=Lattice.from_parameters(*base_entry['parent_metric']);key='dp3_110__'+base_entry['metric_name']+'__1,1,0'
        cache=JointCache(base,repeats,lattice,dataset,ed['registrations'][key],affine['corrections'][key],*raw[mask].T,cfg['width_deg'],cfg['wavelength_A'],nonnegative_background=cfg.get('nonnegative_background',False))
        for weight in cfg['ed_weights']:
            for seed in cfg['seeds']:
                identifier=base_entry['id']+f'__6c_joint_w{weight}_seed{seed}';folder=out/identifier;folder.mkdir(exist_ok=True)
                best,log=search(cache,seed,cfg['steps'],weight,folder/'checkpoint.json',signature);s=cache.to_structure(best)
                angles,intensity=powder(s,cfg['wavelength_A']);np.testing.assert_allclose(angles,cache.angles,atol=1e-7);np.testing.assert_allclose(intensity,cache.intensities(cache.xrd.rebuild(best)),rtol=1e-7,atol=1e-5)
                layer_counts={str(k):np.bincount(best[slots],minlength=3).tolist() for k,slots in enumerate(cache.layer_slots)}
                for slots in cache.layer_slots:np.testing.assert_array_equal(np.bincount(best[slots],minlength=3),np.bincount(cache.labels[slots],minlength=3))
                checks=validate(s,2*np.prod(repeats));s.to(filename=str(folder/'model.cif'),fmt='cif',symprec=None)
                scores={}
                for name,(lo,hi) in {'superstructure':(20.4,34.6),'high_angle':(35,45),'low_angle':(17,20.3)}.items():
                    m=(raw[:,0]>=lo)&(raw[:,0]<=hi)
                    if cfg.get('nonnegative_background',False):
                        v=NonnegativeProfileObjective(angles,*raw[m].T,cfg['width_deg']).evaluate(intensity)
                        v['R2_vs_linear_background']=v['R2_vs_nonnegative_linear_background']
                    else:v=ProfileObjective(angles,*raw[m].T,cfg['width_deg']).evaluate(intensity)
                    scores[name]={k:value for k,value in v.items() if k!='profile'}
                edscores=[]
                for d in ed['datasets']:
                    zone=[1,0,-1] if d['id']=='DP01_101' else d['nominal_zone'];k=d['id']+'__'+base_entry['metric_name']+'__'+','.join(map(str,zone));reg=ed['registrations'][k]
                    hk,hp,g,xy=plane_reflections_3d(lattice,repeats,zone,14);ei=electron_intensity(s,hk,np.linalg.norm(g,axis=1));pp=transform(xy,reg['parameters'],reg['mirror'])
                    write_json(folder/(d['id']+'_reflections.json'),{'hkl_supercell':hk.tolist(),'hkl_parent':hp.tolist(),'g_similarity_nm_inv':pp.tolist(),'kinematic_intensity':ei.tolist()})
                    for mode in ['similarity','affine_diagnostic']:
                        pred=pp if mode=='similarity' else apply_correction(pp,affine['corrections'][k]);visible=in_view(pred,d);maximum=max(float(ei[visible].max()),1e-30)
                        for threshold in [.1,1,5]:
                            selected=visible&(ei>=maximum*threshold/100);m=match_one_to_one(d['observed'],pred[selected],d['position_tolerances_nm_inv'][1])
                            edscores.append({'dataset':d['id'],'mode':mode,'threshold_percent':threshold,**{k:v for k,v in m.items() if k!='pairs'}})
                row={'id':identifier,'nonnegative_background':cfg.get('nonnegative_background',False),'base':base_entry,'repeats':repeats,'ed_objective_weight':weight,'seed':seed,'search_log':log,'layer_species_counts_Li_Sb_Fe_host':layer_counts,'validation':checks,'xrd_scores':scores,'ed_scores':edscores,
                     'final_optimization_scores':cache.score(cache.xrd.rebuild(best),cache.ed.rebuild(best),weight),'geometry_eligible':checks['exact_decimal_composition'] and checks['sixfold_O_coordination_fraction']==1,'gate1':'not established','gate2':'not established','STEM':'not run'}
                write_json(folder/'metadata.json',row);write_json(folder/'xrd.json',{'two_theta_deg':angles.tolist(),'intensity_arbitrary_common_units':intensity.tolist()})
                rows.append({k:v for k,v in row.items() if k not in ['validation','search_log']}|{'geometry_summary':{k:v for k,v in checks.items() if k!='coordination'},'proposal_steps':log['proposal_steps'],'nontrivial_evaluations':log['nontrivial_evaluations']})
                print(identifier,'XRD',round(scores['superstructure']['R2_vs_linear_background'],4),'ED110',[(r['mode'],round(r['f1'],4)) for r in edscores if r['dataset']=='dp3_110' and r['threshold_percent']==1],flush=True)
    summary={'config':cfg,'elapsed_seconds':time.monotonic()-t0,'rows':rows,'validated_shortlist':[],
             'limitations':['This is optimization of exploratory data, not independent validation.','The ED term is a soft position-support objective using simulated intensity saturation at 1% of maximum, not measured ED intensity fitting.','Only [110] DM3 affine-diagnostic geometry enters optimization; all other views and similarity mode are checked afterwards. They were inspected earlier and are not fresh heldouts.','Configured ED weights are exploratory tradeoff settings, not statistical likelihood weights.','Within-layer swaps retain Li:Sb:Fe-host counts and Fe/Ni ratio, with fixed inherited oxygen/interlayer geometry. No energy relaxation.','Finite local searches do not prove optimality; six-c periodicity is a trial and not established by this fit.']}
    write_json(out/'summary.json',summary)
    snapshots=out/'source_snapshot';snapshots.mkdir(exist_ok=True);(snapshots/'joint_occupational.py').write_bytes(Path(__file__).read_bytes());return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
