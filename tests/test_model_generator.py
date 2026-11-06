from pathlib import Path
import numpy as np
import pytest
from pymatgen.core import Structure,Lattice
from lfnso.model_generator import generate,validate


def fixture_parent():
    # Synthetic geometry used for bookkeeping tests, not experimental structure.
    return Structure(Lattice.orthorhombic(5,9,5),['Li','Li',{'Fe':.67,'Ni':.33},{'Fe':.67,'Ni':.33},'Sb','O','Li'],
        [[0,0,0],[.5,.5,0],[0,.2,0],[.5,.7,0],[0,.4,0],[.1,.2,.2],[0,.2,.5]],
        site_properties={'parent_index':list(range(7)),'parent_y':[0,.5,.2,.7,.4,.2,.2],'parent_z':[0,0,0,0,0,.2,.5]})


def test_coupled_model_preserves_composition_and_provenance():
    p=fixture_parent();s=generate(p,ny=4,nz=4,exchange=.25,stacking=True)
    assert s.composition.almost_equals(p.composition*16)
    assert all('parent_z' in t.properties for t in s)
    for t in s:
        if 'Fe' in t.species:assert t.species['Ni']/t.species['Fe']==pytest.approx(.33/.67)


def test_displacement_bound_and_exact_repeat():
    p=fixture_parent();base=generate(p,ny=4);changed=generate(p,ny=4,displacement_A=.02)
    d=base.lattice.get_all_distances(base.frac_coords,changed.frac_coords).diagonal()
    assert d.max()<=.02000001
    assert base.composition==changed.composition


def test_integer_approx_is_explicit_and_reproducible():
    p=fixture_parent();a=generate(p,ny=3,integer_approx=True,seed=4);b=generate(p,ny=3,integer_approx=True,seed=4)
    assert a==b and a.is_ordered
    assert a.composition['Fe']==4 and a.composition['Ni']==2


def test_invalid_exchange_rejected():
    with pytest.raises(ValueError):generate(fixture_parent(),exchange=1.2)


def test_pareto_does_not_rank_numerical_noise():
    from lfnso.model_assessment import nondominated
    rows=[{'id':'a','objectives':[.8,.12,1]},{'id':'b','objectives':[.8-1e-16,.12,1]}]
    assert nondominated(rows)==['a','b']


def test_sb_fe_host_exchange_preserves_exact_decimal_composition():
    # Two Fe hosts and two Sb sites; preserve Fe/Ni co-occupancy on exchanged host fractions.
    p=fixture_parent();p.append('Sb',[.5,.9,0],properties={'parent_index':7,'parent_y':.9,'parent_z':0})
    s=generate(p,ny=7,harmonic=2,host_exchange=1.)
    assert s.composition.almost_equals(p.composition*7)
    assert all(abs(sum(t.species.values())-1)<1e-9 for t in s)
