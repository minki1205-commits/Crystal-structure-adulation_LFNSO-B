"""Bundle the new 3D trials, with explicit legacy-background diagnostic labels."""
import hashlib,json,shutil,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root/'03_structure_model_search';results=root/'results';results.mkdir(exist_ok=True)
files=[]
for folder in ['coupled_search_01','joint_occupational_01','joint_occupational_02','joint_occupational_03','coupled_review_01','coupled_review_02']:
    files.extend((p,p.relative_to(root).as_posix()) for p in (out/folder).rglob('*') if p.is_file())
for pattern in ['src/lfnso/*.py','tests/*.py','configs/coupled*.json','configs/joint_occupational*.json']:
    files.extend((p,p.relative_to(root).as_posix()) for p in root.glob(pattern))
for name in ['README.md','docs/COUPLED_SEARCH.md','docs/ED_MATCHING.md','pyproject.toml','requirements.lock.txt','scripts/package_3d_results.py']:
    files.append((root/name,name))
for relative in ['ed_matching_01/summary.json','ed_metric_diagnostic_01/summary.json',
                 'occupational_search_03/starting_7b_seed17/model.cif','occupational_search_03/conditional_B_7b_seed29/model.cif',
                 'occupational_search_03/conditional_B_11b_seed29/model.cif','occupational_search_03/conditional_B_11b_seed29/xrd.json']:
    p=out/relative;files.append((p,p.relative_to(root).as_posix()))
files.append((out/'raw_xrd_03/LFNSO21-B_profile.csv','inputs/LFNSO21-B_profile.csv'))
extra={}
for filename in ['coupled_search.json','joint_occupational_nonnegative.json','coupled_review_nonnegative.json']:
    c=json.loads((root/'configs'/filename).read_text());c['xrd_csv']='inputs/LFNSO21-B_profile.csv'
    c['output']=c['output'].replace('../03_structure_model_search/','../reproduced/')
    if filename=='coupled_review_nonnegative.json':
        c['staged_summary']=c['staged_summary'].replace('../03_structure_model_search/','../reproduced/')
        c['joint_summaries']=[p.replace('../03_structure_model_search/','../reproduced/') for p in c['joint_summaries']]
    extra['portable_configs/'+filename]=(json.dumps(c,indent=2)+'\n').encode()
extra['REPRODUCE.txt']='''Read 03_structure_model_search/coupled_review_02/report.html first.
Final scoring uses nonnegative structural scale and background endpoints.
Earlier joint_occupational_01/02 and coupled_review_01 permit negative background;
they are preserved diagnostic records, not final physical-background comparisons.
No candidate is a validated structure. No STEM multislice was performed.

From the extracted root with Python 3.12:
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
python -m lfnso.coupled_search portable_configs/coupled_search.json
python -m lfnso.joint_occupational portable_configs/joint_occupational_nonnegative.json
python -m lfnso.coupled_review portable_configs/coupled_review_nonnegative.json

These commands use included XRD numeric data, original candidate CIFs, and previous
ED position/registration summaries. New outputs go to reproduced/.
Raw microscopy images are NOT included; rerunning original ED extraction needs them.
107 new-cycle CIF outputs include 75 staged structures (bases/controls included),
16 initial diagnostic search endpoints and 16 nonnegative-background endpoints.
Three original base CIFs are also supplied as inputs. Sequence classes are not
proven symmetry-distinct atomistic models. Checkpoint/data/source hashes are retained.
'''.encode()
manifest={};archive=results/'LFNSO_3D_joint_search_results.zip';staging=archive.with_suffix('.zip.tmp')
for p,_ in files:
    if not p.is_file():raise FileNotFoundError(p)
with zipfile.ZipFile(staging,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p,name in files:
        content=p.read_bytes();manifest[name]=hashlib.sha256(content).hexdigest();z.writestr(name,content)
    for name,content in extra.items():manifest[name]=hashlib.sha256(content).hexdigest();z.writestr(name,content)
    z.writestr('SHA256SUMS.json',json.dumps(manifest,indent=2))
with zipfile.ZipFile(staging) as z:
    assert z.testzip() is None
    for name,sha in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==sha
    assert len([p for p in z.namelist() if p.endswith('/model.cif')])==110
staging.replace(archive)
shutil.copy2(out/'coupled_review_02/report.html',results/'latest_joint_review.html')
shutil.copy2(out/'coupled_review_02/scores.csv',results/'latest_joint_scores.csv')
records=[]
for p in sorted(results.iterdir()):
    if p.is_file() and p.name!='manifest.json':records.append({'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(results/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'hashed_files':len(manifest),'new_cycle_CIFs':107,'input_CIFs':3},indent=2))
