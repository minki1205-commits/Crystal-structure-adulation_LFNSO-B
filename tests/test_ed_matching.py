"""Synthetic mathematical controls, not fits to experimental observations."""
import numpy as np
import pytest
from pymatgen.core import Structure,Lattice
from pymatgen.analysis.diffraction.tem import TEMCalculator
from lfnso.ed_matching import (plane_reflections,electron_factor,electron_intensity,
                               match_one_to_one,transform,register_parent)
from lfnso.ed_metric_diagnostic import fit_correction,apply_correction


def test_mott_bethe_factors_match_independent_library_atom_calculation():
    s=Structure(Lattice.cubic(5),['Li','Fe','Ni','Sb','O'],[[0,0,0],[.1,0,0],[.2,0,0],[.3,0,0],[.4,0,0]])
    hk=[(1,0,0),(1,1,0),(2,1,1)];d={h:s.lattice.d_hkl(h) for h in hk}
    tem=TEMCalculator();expected=tem.electron_scattering_factors(s,tem.bragg_angles(d))
    for element in ['Li','Fe','Ni','Sb','O']:
        actual=electron_factor(element,10/np.array(list(d.values())))
        np.testing.assert_allclose(actual,[expected[element][h] for h in hk],rtol=1e-12)


def test_partial_occupancies_are_weighted_not_counted_as_full_atoms():
    s=Structure(Lattice.cubic(5),[{'Fe':.67,'Ni':.33}],[[.1,.2,.3]])
    hk=np.array([[1,0,0],[1,1,0]]);g=np.linalg.norm(hk,axis=1)*2
    expected=(.67*electron_factor('Fe',g)+.33*electron_factor('Ni',g))**2
    np.testing.assert_allclose(electron_intensity(s,hk,g),expected)


def test_b_only_supercell_cannot_add_fractional_reflections_in_exact_110_plane():
    lattice=Lattice.monoclinic(5.2,8.9,5.1,109)
    _,hp1,g1,_=plane_reflections(lattice,1,[1,1,0],10)
    _,hp7,g7,_=plane_reflections(lattice,7,[1,1,0],10)
    np.testing.assert_allclose(hp1,hp7);np.testing.assert_allclose(g1,g7)
    _,hp001,_,_=plane_reflections(lattice,7,[0,0,1],10)
    assert np.any(abs(hp001[:,1]-np.rint(hp001[:,1]))>.01)
    for zone in [[1,0,1],[1,0,-1],[0,0,1],[1,1,0]]:
        _,hp,_,_=plane_reflections(lattice,7,zone,10)
        np.testing.assert_allclose(hp@zone,0,atol=1e-12)


def test_exact_repeat_has_same_electron_pattern_and_extinct_fractional_grid():
    p=Structure(Lattice.monoclinic(5.2,8.9,5.1,109),['Li',{'Fe':.67,'Ni':.33},'O'],[[0,0,0],[.5,.2,.3],[.2,.3,.4]])
    a,_,g,_=plane_reflections(p.lattice,1,[0,0,1],10)
    i1=electron_intensity(p,a,np.linalg.norm(g,axis=1))
    s=p.copy();s.make_supercell([1,7,1])
    b,hp,g,_=plane_reflections(p.lattice,7,[0,0,1],10)
    i7=electron_intensity(s,b,np.linalg.norm(g,axis=1))
    integer=np.isclose(hp[:,1],np.rint(hp[:,1]))
    np.testing.assert_allclose(i7[integer]/49,i1,rtol=1e-10,atol=1e-10)
    assert i7[~integer].max()<1e-20


def test_unique_matching_and_cardinality_precede_residual_distance():
    m=match_one_to_one([[0,0],[.12,0]],[[.05,0],[-.08,0]],.1)
    assert m['matches']==2
    m=match_one_to_one([[0,0],[.08,0]],[[.03,0]],.1)
    assert m['matches']==1 and m['recall']==.5
    assert match_one_to_one([],[[0,0]],.1)['matches']==0


def test_dense_grid_does_not_earn_perfect_f1_from_recall_alone():
    observed=np.array([[0,0],[1,0],[2,0]])
    sparse=match_one_to_one(observed,observed,.05)
    dense=match_one_to_one(observed,np.c_[np.arange(0,3,.01),np.zeros(300)],.05)
    assert dense['recall']==1 and sparse['f1']==1 and dense['f1']<.03


def test_parent_similarity_registration_recovers_known_synthetic_mapping():
    lattice=Lattice.monoclinic(5.2,8.9,5.1,109)
    _,hp,_,reference=plane_reflections(lattice,1,[0,0,1],10)
    reference=reference[(hp[:,0]+hp[:,1]).astype(int)%2==0]
    reference=reference[np.argsort(np.linalg.norm(reference,axis=1))]
    observed=transform(reference,[.43,1.03,.012,-.017],-1)
    result=register_parent(observed,reference,[.95,1.05],.04)
    assert result['trimmed_RMS_nm_inv']<1e-5
    assert result['parameters'][1]==pytest.approx(1.03,abs=1e-5)


def test_affine_diagnostic_recovers_small_synthetic_distortion():
    x,y=np.meshgrid(np.arange(-2,3)*1.7,np.arange(-2,3)*2.1)
    reference=np.c_[x.ravel(),y.ravel()]
    matrix=np.array([[1.015,.007],[-.003,.99]])
    observed=reference@matrix.T+[.01,-.02]
    correction=fit_correction(observed,reference,[0,0],.05,.35,.05)
    np.testing.assert_allclose(apply_correction(reference,correction),observed,atol=1e-6)
