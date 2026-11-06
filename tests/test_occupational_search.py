"""Numerical and restart checks on synthetic geometry, not measured structures."""
import json

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from lfnso.occupational_search import DiffractionCache, ProfileObjective, chain


def synthetic_structure():
    return Structure(
        Lattice.monoclinic(5.1, 8.9, 5.2, 109),
        ['Li', 'Li', 'Sb', 'Sb', {'Fe': .67, 'Ni': .33},
         {'Fe': .67, 'Ni': .33}, 'Li', 'O', 'O'],
        [[0, 0, 0], [.5, .5, 0], [0, .2, 0], [.5, .7, 0],
         [0, .4, 0], [.5, .9, 1 - 1e-15], [0, .2, .5],
         [.2, .2, .23], [.8, .8, .77]],
    )


def test_periodic_boundary_slab_and_independent_diffraction():
    cache = DiffractionCache(synthetic_structure())
    assert len(cache.slab) == 6
    assert np.bincount(cache.initial).tolist() == [2, 2, 2]
    check = cache.verify_against_pymatgen(cache.initial)
    assert check['max_relative_intensity_difference'] < 1e-6


def test_incremental_swaps_preserve_composition_and_match_rebuild():
    cache = DiffractionCache(synthetic_structure())
    labels = cache.initial.copy()
    f = cache.factors_for(labels)
    rng = np.random.default_rng(2)
    for _ in range(100):
        i, j = rng.choice(len(labels), 2, replace=False)
        f = cache.swap(f, labels, i, j)
        labels[i], labels[j] = labels[j], labels[i]
    np.testing.assert_allclose(f, cache.factors_for(labels), atol=1e-9)
    assert cache.to_structure(labels).composition == cache.structure.composition
    cache.verify_against_pymatgen(labels)


def test_profile_amplitude_and_background_recovery():
    x = np.linspace(20, 30, 501)
    angles = np.array([22, 25, 28])
    intensity = np.array([2, 3, 1])
    kernel = np.exp(-4 * np.log(2) * ((x[:, None] - angles) / .25)**2) / .25
    y = 2.3 * (kernel @ intensity) + 17 + .4 * (x - x.mean())
    obj = ProfileObjective(angles, x, y, .25)
    result = obj.evaluate(intensity)
    assert result['positive_scale'] == pytest.approx(2.3)
    assert result['linear_background'] == pytest.approx([17, .4])
    assert result['R2_vs_linear_background'] == pytest.approx(1)
    anti = ProfileObjective(angles, x, 40 - y, .25).evaluate(intensity)
    assert anti['positive_scale'] == 0
    assert anti['R2_vs_linear_background'] == pytest.approx(0, abs=1e-10)


def test_interrupted_restart_matches_uninterrupted_chain(tmp_path, monkeypatch):
    import lfnso.occupational_search as module
    cache = DiffractionCache(synthetic_structure())
    x = np.linspace(20, 30, 201)
    y = 100 + 500 * np.exp(-((x - 24.2) / .2)**2)
    obj = ProfileObjective(cache.angles, x, y, .25)
    expected, _ = chain(cache, obj, 12, 500, tmp_path / 'full.json', 'same')
    writer = module.write_json

    def interrupt_after_save(path, data):
        writer(path, data)
        if data['step'] == 250:
            raise InterruptedError('simulated process interruption')

    monkeypatch.setattr(module, 'write_json', interrupt_after_save)
    partial = tmp_path / 'partial.json'
    with pytest.raises(InterruptedError):
        chain(cache, obj, 12, 500, partial, 'same')
    monkeypatch.setattr(module, 'write_json', writer)
    actual, _ = chain(cache, obj, 12, 500, partial, 'same')
    np.testing.assert_array_equal(actual, expected)
    assert json.loads(partial.read_text())['step'] == 500
    again, log = chain(cache, obj, 12, 500, partial, 'same')
    np.testing.assert_array_equal(again, expected)
    assert log['new_diffraction_evaluations_this_invocation'] == 0
    with pytest.raises(ValueError, match='Checkpoint'):
        chain(cache, obj, 12, 500, partial, 'changed')


def test_spectrum_distinguishes_repeat_cell_from_modulation_period():
    from lfnso.occupational_assessment import occupational_spectrum
    parent = synthetic_structure()
    parent.add_site_property('parent_index', list(range(len(parent))))
    parent.add_site_property('parent_y', parent.frac_coords[:, 1].tolist())
    parent.make_supercell([1, 4, 1])
    cache = DiffractionCache(parent)
    labels = cache.initial.copy()
    for j, index in enumerate(cache.slab):
        site = parent[index]
        tile = int(round(site.frac_coords[1] * 4 - site.properties['parent_y'])) % 4
        if tile % 2 and labels[j] < 2:
            labels[j] = 1 - labels[j]
    spectrum = occupational_spectrum(parent, cache, labels, 4)
    assert spectrum[0]['normalized_occupational_power'] == pytest.approx(0)
    assert spectrum[1]['normalized_occupational_power'] == pytest.approx(1)
    assert spectrum[1]['fraction_bstar'] == .5
