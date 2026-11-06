import numpy as np
from lfnso.joint_occupational import JointCache,search
from lfnso.occupational_search import ProfileObjective
from lfnso.coupled_search import powder
from lfnso.coupled_models import plane_reflections_3d
from pymatgen.core import Lattice,Structure
from lfnso.model_generator import FE_NI


def fixture_cache():
    parent=Structure(Lattice.monoclinic(5.2,8.9,5.1,109),['Li','Li','Sb','Sb',FE_NI,FE_NI,'Li','O','O'],
                     [[0,0,0],[.5,.5,0],[0,.2,0],[.5,.7,0],[0,.4,0],[.5,.9,1-1e-15],[0,.2,.5],[.2,.2,.23],[.8,.8,.77]])
    structure=parent.copy();structure.make_supercell([1,1,2]);_,_,_,xy=plane_reflections_3d(parent.lattice,[1,1,2],[1,1,0],10)
    observed=xy[np.argsort(np.linalg.norm(xy,axis=1))][:20]
    dataset={'observed':observed.tolist(),'g_range':[.6,12],'pixel_scale':.1,'center':[500,500],'roi':[0,0,1000,1000]}
    registration={'parameters':[0,1,0,0],'mirror':1};correction={'matrix':np.eye(2).tolist(),'base_offset':[0,0],'total_offset':[0,0]}
    x=np.linspace(20.4,34.6,201);y=100+100*np.exp(-((x-24.5)/.2)**2)
    cache=JointCache(structure,[1,1,2],parent.lattice,dataset,registration,correction,x,y,.25,1.5406)
    return cache,x,y


def test_full_3d_cache_updates_match_rebuilt_structure_and_profile_objective():
    cache,x,y=fixture_cache();labels=cache.labels.copy();fx=cache.xrd.rebuild(labels);fe=cache.ed.rebuild(labels)
    rng=np.random.default_rng(7)
    for _ in range(60):
        slots=cache.layer_slots[rng.integers(2)];i,j=rng.choice(slots,2,replace=False)
        fx=cache.xrd.swap(fx,labels,i,j);fe=cache.ed.swap(fe,labels,i,j);labels[i],labels[j]=labels[j],labels[i]
    np.testing.assert_allclose(fx,cache.xrd.rebuild(labels),atol=1e-9)
    np.testing.assert_allclose(fe,cache.ed.rebuild(labels),atol=1e-9)
    angles,intensity=powder(cache.to_structure(labels))
    np.testing.assert_allclose(angles,cache.angles,atol=1e-7)
    np.testing.assert_allclose(intensity,cache.intensities(fx),atol=1e-6,rtol=1e-8)
    expected=ProfileObjective(angles,x,y,.25).correlation(intensity)
    assert abs(cache.score(fx,fe,.4)[1]-expected)<1e-10


def test_layer_counts_and_finished_checkpoint_resume(tmp_path):
    cache,_,_=fixture_cache();path=tmp_path/'checkpoint.json'
    labels,log=search(cache,9,500,.4,path,'same')
    for slots in cache.layer_slots:
        np.testing.assert_array_equal(np.bincount(labels[slots],minlength=3),np.bincount(cache.labels[slots],minlength=3))
    again,again_log=search(cache,9,500,.4,path,'same')
    np.testing.assert_array_equal(again,labels)
    assert again_log['nontrivial_evaluations']==log['nontrivial_evaluations']
    assert cache.to_structure(labels).composition==cache.structure.composition
