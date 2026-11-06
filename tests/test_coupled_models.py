"""Mathematical invariants and negative controls for c-layer hypotheses."""
import numpy as np
import pytest
from pymatgen.core import Structure,Lattice
from lfnso.model_generator import FE_NI
from lfnso.coupled_models import plane_reflections_3d,primitive_binary_necklaces,extend_layers
from lfnso.coupled_search import powder,verify_powder
from lfnso.ed_matching import plane_reflections


def base_fixture():
    # Bookkeeping-only geometry with two of each slab species and fixed interlayer Li/O.
    return Structure(Lattice.monoclinic(5.2,8.9,5.1,109),['Li','Li','Sb','Sb',FE_NI,FE_NI,'Li','O','O'],
                     [[0,0,0],[.5,.5,0],[0,.2,0],[.5,.7,0],[0,.4,0],[.5,.9,1-1e-15],[0,.2,.5],[.2,.2,.23],[.8,.8,.77]])


def test_general_plane_agrees_with_b_only_and_all_zone_laws():
    lattice=base_fixture().lattice
    for zone in [[0,0,1],[1,0,-1],[1,1,0]]:
        old=plane_reflections(lattice,7,zone,8)
        new=plane_reflections_3d(lattice,[1,7,1],zone,8)
        for a,b in zip(old,new):np.testing.assert_allclose(a,b,atol=1e-12)
        hk,hp,g,xy=plane_reflections_3d(lattice,[1,7,6],zone,8)
        np.testing.assert_allclose(hp@zone,0,atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(g,axis=1),np.linalg.norm(xy,axis=1),atol=1e-12)
    _,hp,_,_=plane_reflections_3d(lattice,[1,7,6],[1,1,0],8)
    assert np.any(abs(hp[:,2]-np.rint(hp[:,2]))>.01)
    np.testing.assert_allclose(hp[:,1],np.rint(hp[:,1]),atol=1e-12)


def test_binary_sequences_are_primitive_and_only_rotation_quotiented():
    assert [len(primitive_binary_necklaces(n)) for n in [2,3,4,6]]==[1,2,3,9]
    assert (0,0,1) in primitive_binary_necklaces(3)
    assert (0,1,1) in primitive_binary_necklaces(3)
    assert (0,1,0,1) not in primitive_binary_necklaces(4)


@pytest.mark.parametrize('sequence,operation', [([0,1],'swap_Sb_Fe'),([0,0,1,1],'swap_Sb_Fe'),([0,1,2],'cycle_Li_Sb_Fe')])
def test_c_layer_operations_conserve_composition_and_fe_ni_ratio(sequence,operation):
    p=base_fixture();s=extend_layers(p,sequence,operation)
    assert s.composition.almost_equals(p.composition*len(sequence))
    for site in s:
        if 'Fe' in site.species:assert site.species['Ni']/site.species['Fe']==pytest.approx(.33/.67)
    # Includes the slab site represented just below fractional z=1.
    assert sum('Fe' in site.species for site in s)==2*len(sequence)


def test_registry_uses_explicit_parent_repeat_and_moves_oxygen_with_slab():
    p=base_fixture();p.make_supercell([1,7,1]);unchanged=extend_layers(p,[0,0,1,1],'none',base_b_repeat=7)
    moved=extend_layers(p,[0,0,1,1],'none',registry_fraction=1/48,base_b_repeat=7)
    d=unchanged.lattice.get_all_distances(unchanged.frac_coords,moved.frac_coords).diagonal()
    assert d.max()==pytest.approx(8.9/48)
    assert moved.composition==unchanged.composition
    for i,site in enumerate(unchanged):
        if site.species_string=='O':assert min(abs(d[i]),abs(d[i]-8.9/48))<1e-8
        if site.species_string=='Li' and abs(site.properties['base_z']-.5)<1e-7:
            assert min(abs(d[i]),abs(d[i]-8.9/96),abs(d[i]-8.9/48))<1e-8


def test_breathing_leaves_oxygen_and_interlayer_li_fixed():
    p=base_fixture();base=extend_layers(p,[0,0,0],'none');moved=extend_layers(p,[0,0,0],'none',breathing_A=.04)
    d=base.lattice.get_all_distances(base.frac_coords,moved.frac_coords).diagonal()
    assert d.max()==pytest.approx(.08)
    for i,site in enumerate(base):
        if site.species_string=='O' or (site.species_string=='Li' and abs(site.properties['base_z']-.5)<1e-7):assert d[i]<1e-8


def test_static_powder_crosscheck_and_exact_repeat_extinctions():
    p=base_fixture();angles,intensity=powder(p);verify_powder(p,angles,intensity,1.5406)
    repeated=extend_layers(p,[0,0],'none');aa,ii=powder(repeated)
    verify_powder(repeated,aa,ii,1.5406)
    m=intensity/intensity.max()*100>1e-5;n=ii/ii.max()*100>1e-5
    np.testing.assert_allclose(angles[m],aa[n],atol=1e-7)
    np.testing.assert_allclose(intensity[m]/intensity.max(),ii[n]/ii.max(),atol=1e-9)
    changed=extend_layers(p,[0,1],'swap_Sb_Fe');a,i=powder(changed);verify_powder(changed,a,i,1.5406)
