import numpy as np
import pytest
from pymatgen.core import Lattice
from lfnso.position_screening import empirical_satellites, xrd_screen


def synthetic_peaks(scale=1):
    # Synthetic reciprocal lattice used only to verify recovery and unit invariance.
    basis=np.array([[3.,.7],[.24,-1.04]])*scale
    hk=[(h,k) for h in range(-2,3) for k in range(-4,5) if (h+k)%2==0 and (h,k)!=(0,0)]
    points=[{'g_yx_nm_inv':(np.array(p)@basis).tolist(),'background_subtracted_peak':200} for p in hk]
    for h,k in hk:
        points.append({'g_yx_nm_inv':(np.array([h,k+.29])@basis).tolist(),'background_subtracted_peak':80})
    return {'peaks':points},[hk.index((1,1)),hk.index((-1,1))]


def test_empirical_q_is_independent_of_detector_scale():
    p,seeds=synthetic_peaks();a=empirical_satellites(p,seeds)
    p,seeds=synthetic_peaks(1.3);b=empirical_satellites(p,seeds)
    assert a['median_delta_v']==pytest.approx(.29)
    assert b['median_delta_v']==pytest.approx(.29)
    assert b['q_norm_nm_inv']==pytest.approx(a['q_norm_nm_inv']*1.3)
    assert a['parent_count']==22


def test_position_screening_retains_fixed_metric_qualification():
    lattice=Lattice.monoclinic(5.1706,8.9382,5.1635,109.49)
    result=xrd_screen(lattice,[[1,1,0],[-1,1,1],[0,2,1]],[4.03,3.62,3.08],[1/6,.25,1/3],.02)
    assert result['discrete_hypotheses'][0]['delta_bstar']==.25
    assert result['discrete_hypotheses'][0]['within_manual_reading_intervals']==3
    assert 'fixed metric' in result['basis_assumptions']


def test_monoclinic_metric_identity_independent_of_a_and_beta():
    for beta in [95,109.49,120]:
        l=Lattice.monoclinic(5.17,8.94,5.16,beta)
        assert l.d_hkl([0,2,1])**-2==pytest.approx(l.d_hkl([0,2,0])**-2+l.d_hkl([0,0,1])**-2)
