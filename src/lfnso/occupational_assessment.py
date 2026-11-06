"""Post-search diagnostics. No post-hoc score constitutes a structure-validation gate."""
import argparse
import base64
import csv
import html
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser

from .model_generator import make_parent
from .occupational_search import DiffractionCache, ProfileObjective
from .pipeline import digest, write_json


def occupational_spectrum(parent, cache, labels, period):
    """Species one-hot Fourier power by parent slot; not a diffraction intensity."""
    parent_ids = sorted({parent[i].properties['parent_index'] for i in cache.slab})
    grid = np.zeros((len(parent_ids), period, 3))
    for i, label in zip(cache.slab, labels):
        site = parent[i]
        tile = int(round(site.frac_coords[1] * period - site.properties['parent_y'])) % period
        grid[parent_ids.index(site.properties['parent_index']), tile, label] = 1
    if not np.all(grid.sum(axis=2) == 1):
        raise ValueError('Incomplete slot-to-tile mapping')
    spectrum = abs(np.fft.fft(grid - grid.mean(axis=1, keepdims=True), axis=1))**2
    power = spectrum.sum(axis=(0, 2))
    folded = np.array([power[m] + (power[-m] if 2*m != period else 0)
                       for m in range(1, period//2 + 1)])
    total = float(folded.sum())
    return [{'m': m, 'fraction_bstar': m/period, 'normalized_occupational_power': float(p/total) if total else 0}
            for m, p in enumerate(folded, 1)]


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    cp = Path(config_path).resolve()
    cfg = json.loads(cp.read_text())
    root = cp.parent
    out = (root / cfg['output']).resolve()
    out.mkdir(parents=True, exist_ok=False)
    data = np.loadtxt(cfg['experimental_csv'], delimiter=',', skiprows=1)
    windows = {'full_superstructure': (20.4, 34.6), 'late_superstructure': (30.3, 34.6),
               'high_angle': (35, 45), 'low_angle_anchors': (17, 20.3), 'common_scale_17_45': (17, 45)}
    rows, curves, sources, logs = [], {}, {}, []
    for source in cfg['sources']:
        folder = (root / source).resolve()
        summary = json.loads((folder / 'summary.json').read_text())
        sources[folder.name] = {'summary_sha256': digest(folder/'summary.json'),
                                'search_provenance': json.loads((folder/'provenance.json').read_text())}
        for row in summary['rows']:
            model_folder = folder / row['id']
            metadata = json.loads((model_folder/'metadata.json').read_text())
            saved = json.loads((model_folder/'xrd.json').read_text())
            intensity = np.array(saved['intensity_arbitrary_common_units'])
            angles = np.array(saved['two_theta_deg'])
            period = row['period_b']
            parent = make_parent(Structure.from_file(summary['config']['template']), [8, 9], row['metric_parameters'])
            parent.make_supercell([1, period, 1])
            cache = DiffractionCache(parent, summary['config']['wavelength_A'])
            np.testing.assert_allclose(angles, cache.angles, atol=1e-8)
            labels = metadata['labels_Li_Sb_FeNi']
            np.testing.assert_allclose(intensity, cache.intensities(cache.factors_for(labels)), rtol=1e-8, atol=1e-7)
            # Ensure the actual exported CIF, not only the in-memory configuration, is checked.
            exported = CifParser(model_folder/'model.cif', frac_tolerance=0).parse_structures(primitive=False)[0]
            exported_cache = DiffractionCache(exported, summary['config']['wavelength_A'])
            exported_cache.verify_against_pymatgen(exported_cache.initial)
            np.testing.assert_allclose(exported_cache.angles, angles, atol=1e-7)
            cif_i = exported_cache.intensities(exported_cache.factors_for(exported_cache.initial))
            np.testing.assert_allclose(cif_i / cif_i.max(), intensity / intensity.max(), atol=2e-6)
            record = {'cycle': folder.name, 'id': row['id'], 'period_b': period,
                      'metric': row['metric_name'], 'site_count': len(parent),
                      'exact_composition': row['geometry_summary']['exact_decimal_composition'],
                      'CN6_fraction': row['geometry_summary']['sixfold_O_coordination_fraction'],
                      'cif_roundtrip_diffraction_pass': True, 'cif_reader_fraction_rationalization_tolerance': 0,
                      'occupational_spectrum': occupational_spectrum(parent, cache, labels, period), 'width_scores': {}}
            for width in cfg['widths_deg']:
                per_window, per_control = {}, {}
                for name, (lo, hi) in windows.items():
                    mask = (data[:, 0] >= lo) & (data[:, 0] <= hi)
                    obj = ProfileObjective(angles, *data[mask].T, width)
                    fit = obj.evaluate(intensity)
                    baseline = obj.evaluate(cache.intensities(cache.factors_for(cache.initial)))
                    per_window[name] = {k:v for k,v in fit.items() if k != 'profile'}
                    per_control[name] = {k:v for k,v in baseline.items() if k != 'profile'}
                    if width == .25:
                        curves[(folder.name, row['id'], name)] = (obj.x, obj.y, fit['profile'], baseline['profile'])
                record['width_scores'][str(width)] = {'candidate': per_window, 'ordered_control': per_control}
            rows.append(record)
            logs.append(metadata['search_log'])
    # Descriptive ordering at the pre-specified search width only; no validated shortlist.
    target_cycle = (root / cfg['sources'][-1]).resolve().name
    ranked = sorted([r for r in rows if r['cycle'] == target_cycle],
                    key=lambda r:r['width_scores']['0.25']['candidate']['full_superstructure']['R2_vs_linear_background'], reverse=True)
    summary = {'config': cfg, 'sources': sources,
               'assessment_code_sha256': digest(Path(__file__)),
               'dependencies_sha256': {p.name:digest(p) for p in [Path(__file__).with_name('occupational_search.py'), Path(__file__).with_name('model_generator.py')]},
               'raw_profile_sha256': digest(Path(cfg['experimental_csv'])),
               'proposal_steps_total': sum(r['proposal_steps'] for r in logs),
               'diffraction_evaluations_recorded': sum(r['new_diffraction_evaluations_this_invocation'] for r in logs),
               'candidate_count': len(rows), 'rows': rows,
               'descriptive_best_full_superstructure_at_0p25': ranked[0]['id'],
               'validated_shortlist': [], 'gate_status': 'not established; visible unexplained/extra peaks, no calibrated ED validation',
               'interpretation': ['All are unrelaxed trial structures, exact average Fe/Ni composition, fixed interlayer Li and O.',
                                  'Width sensitivity is post-search evaluation of fixed labels, not re-optimized structures.',
                                  'Each window has its own nonnegative scale and linear background. R2 is improvement over a linear-background-only fit, not Rwp or a probability.',
                                  '35-45 degree data and the failed 30.3-34.6 degree window informed follow-up design; no fresh independent validation is claimed.',
                                  'Supercell periodicity allows many harmonics; occupational Fourier power is not electron diffraction intensity.',
                                  'Candidate improvement at high angle must be compared with its ordered control under the same metric.',
                                  'Oxygen/cation relaxation, stacking/disorder, instrument profile and full-pattern refinement remain unresolved.']}
    write_json(out/'summary.json', summary)
    snapshots = out/'source_snapshot'
    snapshots.mkdir()
    for source in [Path(__file__), Path(__file__).with_name('occupational_search.py'), Path(__file__).with_name('model_generator.py'), cp]:
        (snapshots/source.name).write_bytes(source.read_bytes())
    with (out/'scores.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['cycle','id','width_deg','window','candidate_R2_vs_linear_background','control_R2_vs_linear_background','candidate_normalized_RMS'])
        for r in rows:
            for width, scores in r['width_scores'].items():
                for window, score in scores['candidate'].items():
                    w.writerow([r['cycle'],r['id'],width,window,score['R2_vs_linear_background'],scores['ordered_control'][window]['R2_vs_linear_background'],score['normalized_RMS_by_data_range']])
    best = ranked[0]
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    for ax, window in zip(axes.flat, ['full_superstructure','late_superstructure','high_angle','low_angle_anchors']):
        x, y, fit, baseline = curves[(best['cycle'],best['id'],window)]
        ax.plot(x,y,'k-',lw=.8,label='LFNSO21-B measured')
        ax.plot(x,fit,'r-',lw=1,label='Trial arrangement + linear baseline')
        ax.plot(x,baseline,color='steelblue',lw=.8,alpha=.8,label='Ordered control + linear baseline')
        scores = best['width_scores']['0.25']
        ax.set(title=f"{window}: trial R2={scores['candidate'][window]['R2_vs_linear_background']:.3f}, control={scores['ordered_control'][window]['R2_vs_linear_background']:.3f}",xlabel='2theta (degrees)',ylabel='Intensity (input units)')
        ax.legend(fontsize=7)
    fig.suptitle(best['id']+'; diagnostic comparison, separate scale/background per window')
    fig.tight_layout(); fig.savefig(out/'best_diagnostic.png',dpi=150); plt.close(fig)
    notes = (f'Li₄Fe₀.₆₇Ni₀.₃₃SbO₆ 조성 및 Fe/Ni 공동 점유를 유지한 {len(rows)}개 탐색 종점입니다. '
             f"최종 구조가 아닙니다. {summary['proposal_steps_total']:,}회 swap proposal에는 같은 원소 선택으로 건너뛴 단계와 반복 방문이 포함되며 서로 다른 구조 개수가 아닙니다. "
             '4b/7b/11b 및 두 격자를 시험했습니다. 아래 R²는 선형 배경만 맞춘 경우 대비 개선량입니다. '
             '높은 각도에서는 조건부 격자 자체의 개선이 크므로 ordering의 증거로 단독 해석할 수 없습니다. '
             '34.35° 피크를 포함하여 다시 탐색하면 해당 피크를 만들 수 있지만 불필요한 반사와 강도 불일치가 남습니다. '
             'Cu Kα 1.5406 Å는 잠정값이고, 산소 좌표·격자는 최종 정련하지 않았습니다.')
    page = '<!doctype html><meta charset="utf-8"><title>LFNSO21-B occupancy search assessment</title><style>body{font:16px system-ui;max-width:1250px;margin:30px auto}img{max-width:100%}td,th{padding:7px;border-bottom:1px solid #ddd}pre{white-space:pre-wrap}</style>'
    page += '<h1>자리 배열 탐색 결과 · 검증 전 후보</h1><p>'+notes+'</p>'
    page += '<img src="data:image/png;base64,'+base64.b64encode((out/'best_diagnostic.png').read_bytes()).decode()+'">'
    page += f'<h2>추가 탐색 {len(ranked)}개 후보 · 폭 0.25°</h2><table><tr><th>후보</th><th>20.4–34.6° R²</th><th>35–45° R²</th><th>35–45° ordered control R²</th></tr>'
    for r in ranked:
        s = r['width_scores']['0.25']
        page += '<tr><td>'+r['id']+'</td>'+''.join(f'<td>{v:.4f}</td>' for v in [s['candidate']['full_superstructure']['R2_vs_linear_background'],s['candidate']['high_angle']['R2_vs_linear_background'],s['ordered_control']['high_angle']['R2_vs_linear_background']])+'</tr>'
    page += '</table><p>16–46° XRD 계산을 pymatgen으로 교차 확인했고, 내보낸 CIF를 다시 읽어 회절을 검증했습니다. CN6는 2.7 Å cutoff 아래 기하 검사이며 에너지 안정성 검증이 아닙니다. 전체 점수, 피크 폭 민감도, Fourier 성분과 provenance는 아래에 있습니다.</p><pre>'+html.escape(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre>'
    (out/'report.html').write_text(page)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    run(parser.parse_args().config)
