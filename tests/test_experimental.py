import numpy as np
import pytest
from lfnso.experimental import calibrated_fft,extract_peaks


def test_fft_frequency_calibration():
    n=128;pixel=.02
    y,x=np.indices((n,n));a=np.cos(2*np.pi*10*x/n)
    f,dq=calibrated_fft(a,[pixel,pixel],window=False)
    assert dq[1]==pytest.approx(1/(n*pixel))
    peaks=np.argpartition(np.abs(f).ravel(),-2)[-2:]
    coords=np.array(np.unravel_index(peaks,f.shape)).T
    assert sorted(coords[:,1].tolist())==[n//2-10,n//2+10]
    assert np.max(np.abs(f.imag))<1e-10


def test_peak_coordinates_center_and_calibration():
    y,x=np.indices((256,256));cy,cx=116.3,138.2
    a=np.zeros_like(y,dtype=float)+5
    for dy,dx,amp in [(0,0,1000),(30,20,200),(-30,-20,200),(20,-40,180),(-20,40,180),(50,30,100),(-50,-30,100)]:
        a+=amp*np.exp(-((y-cy-dy)**2+(x-cx-dx)**2)/(2*2**2))
    r=extract_peaks(a,[.02,.03])
    assert r['opposite_pairs']>=3
    assert np.allclose(r['center_yx_px'],[cy,cx],atol=.15)
    points=np.array([q['g_yx_nm_inv'] for q in r['peaks']])
    assert np.min(np.linalg.norm(points-[.6,.6],axis=1))<.01


def test_nonfinite_input_is_rejected_not_silently_repaired():
    a=np.ones((64,64));a[0,0]=np.nan
    with pytest.raises(ValueError):extract_peaks(a,[.1,.1])
    with pytest.raises(ValueError):calibrated_fft(a,[.1,.1])
