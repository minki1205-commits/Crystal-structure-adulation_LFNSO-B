"""Small executable audit → reference controls → XRD → honest availability report."""
import argparse
import hashlib
import html
import json
from pathlib import Path
import platform
import time
import numpy as np
from pymatgen.core import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from .geometry import reciprocal_cartesian, zone_law


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inventory(roots):
    rows = []
    for name in roots:
        root = Path(name).resolve()
        if not root.exists():
            rows.append({'path': str(root), 'status': 'unavailable'})
            continue
        files = [root] if root.is_file() else sorted(p for p in root.rglob('*') if p.is_file() and '.git' not in p.parts)
        for p in files:
            rows.append({'path': str(p), 'format': p.suffix, 'bytes': p.stat().st_size,
                         'sha256': digest(p), 'sample': None, 'zone_axis': None,
                         'calibration': None, 'experimental_conditions': None,
                         'status': 'unreviewed'})
    return rows


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def run(config_path):
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text())
    # Paths are relative to configuration file, independent of caller cwd.
    resolve = lambda p: (config_path.parent / p).resolve()
    out = resolve(cfg['output'])
    # Never overwrite a previous cycle or an input tree.
    for name in cfg.get('input_roots', []):
        root = resolve(name)
        if out == root or out.is_relative_to(root) or root.is_relative_to(out):
            raise ValueError('Output must be separate from input roots')
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    for directory in ['inputs_manifest','reference_structures','candidate_cifs','xrd_results','ed_fft_results','rankings','stem_shortlist','reports','configs','logs']:
        (out / directory).mkdir()
    write_json(out/'configs/run.json', cfg)
    roots = [str(resolve(p)) for p in cfg.get('input_roots', [])]
    manifest = inventory(roots)
    write_json(out/'inputs_manifest/files.json', manifest)
    report = {'status': 'blocked_missing_parent', 'experimental_constraints': {
        'xrd': 'unavailable', 'ed': 'unavailable', 'fft': 'unavailable', 'stem': 'unavailable'},
        'ranked_candidates': [], 'candidates': [], 'gate1': 'not_evaluated',
        'gate2': 'not_evaluated', 'stem': 'not_run', 'python': platform.python_version(),
        'limitations': ['Experimental scoring, occupational exchange, slab-aware stacking, disorder ensembles and multislice are not implemented in this first increment.']}
    parent_name = cfg.get('parent_cif')
    if parent_name:
        parent_path = resolve(parent_name)
        parent = Structure.from_file(parent_path)
        parent_hash = digest(parent_path)
        report['parent'] = {'path': str(parent_path), 'sha256': parent_hash,
                            'composition': parent.composition.as_dict(), 'lattice': parent.lattice.as_dict(),
                            'role': cfg.get('parent_role', 'user-supplied structure; suitability unverified')}
        # Starting controls only. No inferred site assignments or chemical valences.
        for period in cfg.get('periods', [1,3,4]):
            if not isinstance(period, int) or period < 1 or period > 12:
                raise ValueError('period must be integer from 1 to 12')
            candidate = parent.copy()
            candidate.make_supercell([1, period, 1])
            name = f'exact_repeat_{period}b'
            candidate.to(filename=str(out/'candidate_cifs'/f'{name}.cif'), fmt='cif')
            distances = candidate.distance_matrix.copy()
            np.fill_diagonal(distances, np.inf)
            minimum = float(distances.min()) if len(candidate)>1 else None
            valid_occupancy = all(0 < sum(site.species.values()) <= 1+1e-8 for site in candidate)
            row = {'id': name, 'kind': 'negative_control', 'parent_sha256': parent_hash,
                   'transformation': np.diag([1,period,1]).tolist(), 'seed': cfg.get('seed',0),
                   'composition': candidate.composition.as_dict(), 'occupancy_valid': valid_occupancy,
                   'minimum_distinct_site_distance_A': minimum,
                   'chemical_validity': 'incomplete_requires_valence_and_coordination',
                   'physical_modulation': False, 'gate1': 'not_evaluated', 'gate2': 'not_evaluated'}
            wavelength = cfg.get('wavelength_A')
            if wavelength is not None:
                if not np.isfinite(wavelength) or wavelength <= 0:
                    raise ValueError('Explicit positive wavelength in Å required')
                pattern = XRDCalculator(wavelength=wavelength).get_pattern(candidate, two_theta_range=tuple(cfg.get('two_theta_range',[5,80])))
                write_json(out/'xrd_results'/f'{name}.json', {'two_theta_deg':pattern.x.tolist(), 'relative_intensity':pattern.y.tolist(), 'wavelength_A': wavelength, 'interpretation':'ideal periodic powder prediction; not an experimental fit'})
            reflections = []
            for h in range(-2,3):
                for k in range(-2*period,2*period+1):
                    for l in range(-2,3):
                        if (h,k,l)==(0,0,0): continue
                        parent_hkl = [h,k/period,l]
                        reflections.append({'super_hkl':[h,k,l], 'parent_hkl':parent_hkl,
                          'g_A_inv':reciprocal_cartesian([h,k,l],candidate.lattice.matrix).tolist(),
                          'exact_repeat_translation_allowed': k % period == 0,
                          'parent_zones':[u for u in [[0,0,1],[1,0,1],[1,1,0]] if zone_law(parent_hkl,u)]})
            write_json(out/'ed_fft_results'/f'{name}_geometry.json', {'interpretation':'geometric reflection table, no electron scattering intensities', 'reflections': reflections})
            report['candidates'].append(row)
        report['status'] = 'controls_generated_unranked'
    report['elapsed_seconds'] = time.monotonic()-start
    write_json(out/'reports/result.json',report)
    (out/'reports/report.html').write_text('<!doctype html><meta charset="utf-8"><title>LFNSO research status</title><h1>LFNSO screening status</h1><p>No structure has been established. Missing evidence does not pass screening gates.</p><pre>'+html.escape(json.dumps(report,indent=2,ensure_ascii=False))+'</pre>')
    print(json.dumps({'output':str(out),'status':report['status'],'candidates':len(report['candidates']),'ranked':0}))
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('config', help='JSON configuration; relative paths resolve beside this file')
    args=parser.parse_args()
    run(args.config)

if __name__ == '__main__': main()
