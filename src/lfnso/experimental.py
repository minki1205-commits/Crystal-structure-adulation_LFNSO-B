"""Conservative preliminary peak measurement; no automatic hkl or phase assignment."""
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.spatial import cKDTree


def calibrated_fft(image, pixel_nm, window=True):
    image = np.asarray(image, float)
    if image.ndim != 2 or not np.isfinite(image).all(): raise ValueError('Finite 2D image required')
    pixel_nm = np.asarray(pixel_nm,float)
    if pixel_nm.shape != (2,) or not np.isfinite(pixel_nm).all() or (pixel_nm <= 0).any():
        raise ValueError('Positive y,x calibration in nm required')
    a = image-image.mean()
    if window: a = a*np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))
    return np.fft.fftshift(np.fft.fft2(a)), 1/(np.array(a.shape)*pixel_nm)


def extract_peaks(image, scale_yx, is_fft=False, max_peaks=200):
    """Background-subtracted local maxima. Widths are moments, not fitted FWHM/errors."""
    image=np.abs(np.asarray(image)).astype(float)
    scale=np.asarray(scale_yx,float)
    if image.ndim!=2 or not np.isfinite(image).all():raise ValueError('Finite 2D image required')
    if scale.shape!=(2,) or (scale<=0).any() or not np.isfinite(scale).all():raise ValueError('Invalid reciprocal calibration')
    work=np.log1p(image) if is_fft else image
    smooth=gaussian_filter(work,1.)
    residual=smooth-gaussian_filter(work,12.)
    noise=1.4826*np.median(np.abs(residual-np.median(residual)))
    threshold=max(8*noise, float(residual.max())*(0.01 if is_fft else 0.0005))
    maxima=(residual==maximum_filter(residual,size=9)) & (residual>threshold)
    maxima[:10,:]=False;maxima[-10:,:]=False;maxima[:,:10]=False;maxima[:,-10:]=False
    yy,xx=np.where(maxima)
    order=np.argsort(residual[yy,xx])[::-1][:max_peaks]
    points=[]
    for j in order:
        y,x=int(yy[j]),int(xx[j]);rad=5
        patch=np.maximum(residual[y-rad:y+rad+1,x-rad:x+rad+1],0)
        gy,gx=np.mgrid[y-rad:y+rad+1,x-rad:x+rad+1]
        total=patch.sum()
        cy,cx=(patch*gy).sum()/total,(patch*gx).sum()/total
        sy,sx=np.sqrt((patch*(gy-cy)**2).sum()/total),np.sqrt((patch*(gx-cx)**2).sum()/total)
        points.append({'y_px':float(cy),'x_px':float(cx),'background_subtracted_peak':float(residual[y,x]),
                       'moment_sigma_yx_px':[float(sy),float(sx)]})
    if not points:raise ValueError('No peaks detected')
    coords=np.array([[p['y_px'],p['x_px']] for p in points])
    image_center=np.array(image.shape)//2
    if is_fft:
        center=image_center.astype(float);center_source='DFT index zero';pair_count=None;scatter=None
    else:
        central=np.linalg.norm(coords-image_center,axis=1)<min(image.shape)*.2
        if not central.any():raise ValueError('No central-beam candidate')
        central_indices=np.where(central)[0]
        center=coords[central_indices[0]].copy() # descending peak strength, within central search region
        # Geometric opposite-spot pairing; intensity symmetry is not assumed.
        tree=cKDTree(coords)
        pair_count=0;scatter=None
        for _ in range(3):
            distances,indices=tree.query(2*center-coords)
            pairs={(min(i,int(j)),max(i,int(j))) for i,j in enumerate(indices)
                   if i!=j and distances[i]<6 and np.linalg.norm(coords[i]-center)>15}
            if len(pairs)<3:break
            mids=np.array([(coords[i]+coords[j])/2 for i,j in pairs])
            center=np.median(mids,axis=0);pair_count=len(pairs)
            scatter=np.sqrt(np.mean((mids-center)**2,axis=0)).tolist()
        center_source='opposite-spot midpoints' if pair_count>=3 else 'bright central peak; unvalidated'
    for p in points:
        g=(np.array([p['y_px'],p['x_px']])-center)*scale
        p.update({'g_yx_nm_inv':g.tolist(),'g_norm_nm_inv':float(np.linalg.norm(g)),
                  'classification':'unindexed', 'position_uncertainty':'not_estimated'})
    return {'center_yx_px':center.tolist(),'center_method':center_source,'opposite_pairs':pair_count,
            'midpoint_rms_yx_px':scatter, 'threshold':threshold,'noise_estimator':float(noise),
            'coordinate_convention':'array y increases downward, x right; orientation to crystal not fitted',
            'parameters':{'gaussian_sigma_px':1,'background_sigma_px':12,'nms_size_px':9,'max_peaks':max_peaks},
            'limitations':['Preliminary detections can include artifacts and miss weak satellites.',
                          'Moment widths are window-limited; no instrumental broadening correction.',
                          'Positions require calibration/orientation validation before structural scoring.'],
            'peaks':points}
