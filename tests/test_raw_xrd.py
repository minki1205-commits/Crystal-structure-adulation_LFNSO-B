import numpy as np
import pytest
from openpyxl import Workbook
from lfnso.raw_xrd import read_workbook,fit_peak,delta_from_001_020_satellite,sine_coordinate
from lfnso.peak_groups import fit_group


def test_workbook_columns_and_rows_are_preserved(tmp_path):
    path=tmp_path/'input.xlsx';w=Workbook();s=w.active;s.title='data'
    s.append(['A',None,'B',None]);s.append(['2θ (degree)','Intensity','2θ (degree)','Intensity'])
    for angle in [10,10.02,10.04]:s.append([angle,101,angle,202])
    w.save(path)
    values,meta=read_workbook(path,'data')
    assert values['A'].shape==(3,2)
    assert values['B'][:,1].tolist()==[202,202,202]
    assert meta['samples'][1]['excel_columns_1_based']==[3,4]
    assert meta['samples'][0]['rows_1_based']==[3,5]


def test_incomplete_spreadsheet_is_rejected(tmp_path):
    path=tmp_path/'input.xlsx';w=Workbook();s=w.active;s.title='data'
    s.append(['A',None]);s.append(['2θ (degree)','Intensity']);s.append([10,None]);w.save(path)
    with pytest.raises(ValueError,match='Missing data'):read_workbook(path,'data')


def test_peak_center_recovery_with_slope_and_noise():
    rng=np.random.default_rng(4);x=np.arange(21.5,22.5,.02)
    y=200+20*(x-22)+1200*np.exp(-.5*((x-22.137)/.065)**2)+rng.normal(0,4,len(x))
    r=fit_peak(x,y,22.14)
    assert abs(r['center_two_theta_deg']-22.137)<.002
    assert abs(r['observed_envelope_fwhm_deg']-2.35482*.065)<.006


def test_delta_does_not_require_wavelength():
    # Construct three spacings for an orthogonal b-vs-c reciprocal subspace.
    b,c,delta=9.,4.9,.28
    d=[c,b/2,1/np.sqrt(c**-2+(2+delta)**2/b**2)]
    for wavelength in [1.,1.5406,1.79]:
        angles=2*np.rad2deg(np.arcsin(wavelength/(2*np.array(d))))
        assert delta_from_001_020_satellite(angles)==pytest.approx(delta)


def test_empirical_doublet_recovery():
    rng=np.random.default_rng(1);x=np.arange(34.8,35.8,.01)
    y=150+1000*np.exp(-.5*((x-35.13)/.05)**2)+600*np.exp(-.5*((x-35.39)/.07)**2)+rng.normal(0,3,len(x))
    r=fit_group(x,y)
    centers=[p['center_deg'] for p in r['models'][1]['components']]
    assert np.allclose(centers,[35.13,35.39],atol=.003)
    assert r['aicc_two_minus_one'] < -100
