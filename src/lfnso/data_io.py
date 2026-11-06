"""Read DM image data without altering raw files or inventing calibration."""
from pathlib import Path
import numpy as np
from ncempy.io import dm
from .pipeline import digest


def read_dm(path):
    path = Path(path)
    with dm.fileDM(path) as f:
        dataset = f.getDataset(0)  # NCEMpy skips the thumbnail.
        selected = {}
        for key, value in f.allTags.items():
            if any(key.endswith(s) for s in ['Acquisition Date', 'Acquisition Time', 'Signal Name',
                    '.Voltage', '.STEM Camera Length', '.Indicated Magnification', '.exposure (s)']):
                if isinstance(value, np.generic): value = value.item()
                if isinstance(value, (str, int, float)): selected[key] = value
        layout = {'dataset_count_including_thumbnail': int(f.numObjects),
                  'data_type_codes': [int(v) for v in f.dataType]}
    a = dataset['data']
    if a.ndim != 2: raise ValueError('Only 2D images are currently supported')
    row = {'path': str(path.resolve()), 'sha256': digest(path), 'bytes': path.stat().st_size,
           'shape_yx': list(a.shape), 'dtype': str(a.dtype),
           'pixel_size_yx': [float(v) for v in dataset['pixelSize']],
           'pixel_unit_yx': list(dataset['pixelUnit']),
           'dm_origin_yx': [float(v) for v in dataset['pixelOrigin']],
           'acquisition_metadata': selected, **layout}
    finite = np.isfinite(a)
    magnitudes = np.abs(a[finite]).astype(float)
    row['nonfinite_pixels'] = int((~finite).sum())
    row['finite_magnitude_percentiles'] = np.percentile(magnitudes,[0,1,50,99,100]).tolist()
    # This detects catastrophic numeric values, not ordinary diffraction dynamic range.
    limit = max(float(np.percentile(magnitudes, 99))*1e12, 1e20)
    extreme = finite & (np.abs(a) > limit)
    row['catastrophic_value_threshold'] = limit
    row['catastrophic_pixels'] = int(extreme.sum())
    row['quantitative_read_status'] = 'quarantined_numeric_anomalies' if (~finite).any() or extreme.any() else 'readable'
    return a, row
