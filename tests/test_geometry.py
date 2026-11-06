import numpy as np
import pytest
from lfnso.geometry import reflection_to_super,zone_to_super,reciprocal_cartesian,structure_factors,zone_law

@pytest.mark.parametrize('n',[3,4,6,12])
def test_indexing_and_invariance(n):
    s=np.diag([1,n,1])
    a=np.array([[5.17,0,0],[0,8.94,0],[-1.72,0,4.86]])
    assert np.allclose(reflection_to_super([0,2,0],s),[0,2*n,0])
    assert np.allclose(reflection_to_super([1,1,0],s),[1,n,0])
    assert np.allclose(zone_to_super([1,1,0],s)*n,[n,1,0])
    h=[1,2,3]
    assert np.allclose(reciprocal_cartesian(h,a),reciprocal_cartesian(reflection_to_super(h,s),s@a))

@pytest.mark.parametrize('n',[3,4])
def test_repeat_has_no_satellites(n):
    # Synthetic mathematical fixture, never experimental LFNSO data.
    frac=np.array([[.13,(.27+j)/n,.41] for j in range(n)])
    f=structure_factors(frac,[[0,k,0] for k in range(1,n+1)],np.ones(n))
    assert np.allclose(f[:-1],0,atol=1e-12)
    assert abs(f[-1])==pytest.approx(n)

def test_pure_b_modulation_not_in_110_zolz():
    assert not zone_law([0,.25,0],[1,1,0])
    assert zone_law([0,.25,0],[1,0,1])
