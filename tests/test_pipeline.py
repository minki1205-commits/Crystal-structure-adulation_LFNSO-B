import json
import pytest
from pymatgen.core import Lattice,Structure
from lfnso.pipeline import run

def test_missing_data_never_ranks(tmp_path):
    p=tmp_path/'config.json'
    p.write_text(json.dumps({'input_roots':['missing'],'output':'result'}))
    r=run(p)
    assert r['status']=='blocked_missing_parent'
    assert r['ranked_candidates']==[]
    assert r['gate1']=='not_evaluated'
    with pytest.raises(FileExistsError):run(p)

def test_real_cif_roundtrip_and_xrd_controls(tmp_path):
    # Synthetic generic fixture, not an LFNSO structure hypothesis.
    parent=Structure(Lattice.cubic(4),['Na','Cl'],[[0,0,0],[.5,.5,.5]])
    parent.to(filename=str(tmp_path/'test.cif'))
    p=tmp_path/'config.json'
    p.write_text(json.dumps({'parent_cif':'test.cif','output':'result','periods':[1,3,4],'wavelength_A':1.5406}))
    r=run(p)
    assert len(r['candidates'])==3
    assert r['ranked_candidates']==[]
    for n in [1,3,4]:
        s=Structure.from_file(tmp_path/f'result/candidate_cifs/exact_repeat_{n}b.cif')
        assert s.composition==parent.composition*n
    patterns=[json.loads((tmp_path/f'result/xrd_results/exact_repeat_{n}b.json').read_text()) for n in [1,3,4]]
    import numpy as np
    for pat in patterns[1:]:
        assert np.allclose(pat['two_theta_deg'],patterns[0]['two_theta_deg'])
        assert np.allclose(pat['relative_intensity'],patterns[0]['relative_intensity'])
