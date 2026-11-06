import numpy as np
import pytest
from lfnso.physical_profile import NonnegativeBaseline,NonnegativeProfileObjective
from lfnso.occupational_search import ProfileObjective


def test_positive_background_and_scale_are_recovered():
    x=np.linspace(0,10,501);profile=4+np.exp(-((x-5)/.5)**2)
    y=2.3*profile+17+.4*x
    result=NonnegativeBaseline(x,y).fit(profile,details=True)
    assert result['positive_scale']==pytest.approx(2.3)
    assert result['background_endpoint_values']==pytest.approx([17,21])
    assert result['R2_vs_nonnegative_linear_background']==pytest.approx(1)


def test_negative_background_cannot_be_hidden_by_structure_offset():
    x=np.linspace(0,10,501);profile=4+np.exp(-((x-5)/.5)**2)
    y=2.3*profile-4+.1*x
    result=NonnegativeBaseline(x,y).fit(profile,details=True)
    assert min(result['background_endpoint_values'])>=0
    assert result['R2_vs_nonnegative_linear_background']<.999
    assert min(result['profile'])>=0


def test_constraint_agrees_with_legacy_fit_when_background_positive():
    x=np.linspace(20,30,501);angles=np.array([22.,25.,28.]);intensity=np.array([2.,3.,1.])
    kernel=np.exp(-4*np.log(2)*((x[:,None]-angles)/.25)**2)/.25
    y=2.3*(kernel@intensity)+17+.4*(x-25)
    legacy=ProfileObjective(angles,x,y,.25).evaluate(intensity)
    new=NonnegativeProfileObjective(angles,x,y,.25).evaluate(intensity)
    np.testing.assert_allclose(new['profile'],legacy['profile'],atol=1e-9)
