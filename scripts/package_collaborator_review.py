"""Package the collaborator audit and its exact numeric/image inputs."""
import hashlib,json,shutil,zipfile
from pathlib import Path

root=Path(__file__).resolve().parents[1]
cfg=json.loads((root/'configs/collaborator_review.json').read_text())
cp=root/'configs';out=(cp/cfg['output']).resolve();results=root/'results'
summary=json.loads((out/'summary.json').read_text())
ed_path=(cp/cfg['ed_summary']).resolve();ed=json.loads(ed_path.read_text())
files={p.relative_to(root).as_posix():p for p in out.rglob('*') if p.is_file()}
for pattern in ['src/lfnso/*.py','tests/*.py']:
    files.update({p.relative_to(root).as_posix():p for p in root.glob(pattern)})
for name in ['configs/collaborator_review.json','docs/COLLABORATOR_CIF_REVIEW.md','pyproject.toml','requirements.lock.txt','scripts/package_collaborator_review.py']:
    files[name]=root/name
files['inputs/original_collaborator.cif']=(cp/cfg['cif']).resolve()
files['inputs/LFNSO21-B_profile.csv']=(cp/cfg['xrd_csv']).resolve()
files['inputs/xrd_fits.json']=(cp/cfg['xrd_fits']).resolve()
files['inputs/ed_summary.json']=ed_path
portable=dict(cfg,output='../reproduced/collaborator_review',cif='../inputs/original_collaborator.cif',
              xrd_csv='../inputs/LFNSO21-B_profile.csv',xrd_fits='../inputs/xrd_fits.json',ed_summary='../inputs/ed_summary.json',raw_overrides={})
for dataset in ed['datasets']:
    raw=Path(dataset['raw']);name='inputs/microscopy/'+raw.name;files[name]=raw
    portable['raw_overrides'][dataset['id']]='../'+name
    files['inputs/peak_sources/'+dataset['id']+'.json']=Path(dataset['peak_source'])
extra={'portable_configs/collaborator_review.json':(json.dumps(portable,indent=2)+'\n').encode(),
       'REPRODUCE.txt':'''Read docs/COLLABORATOR_CIF_REVIEW.md and the HTML report first.
The supplied CIF is preserved byte-for-byte in inputs/original_collaborator.cif.
Calculation copies and conditional metric fits are distinct output files.
No validated structure or STEM multislice match is claimed.

From this extracted directory with Python 3.12:
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
python -m lfnso.collaborator_review portable_configs/collaborator_review.json

Outputs go to reproduced/collaborator_review; choose a new path for another run.
Two DM3 and two TIFF files are included for audit overlays, along with prior
detected ED positions and numeric B XRD/profile-fit outputs. This is not a backup
of all microscopy/STEM/XRD original files. Earlier peak extraction and profile-fit
pipelines are not rerun by the audit command. Original input and source hashes
are retained; rerun summary config/path hashes may differ with portable paths.
'''.encode()}
manifest={};archive=results/'LFNSO_collaborator_CIF_review.zip';temp=archive.with_suffix('.zip.tmp')
with zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for name,path in sorted(files.items()):
        blob=path.read_bytes();manifest[name]=hashlib.sha256(blob).hexdigest();z.writestr(name,blob)
    for name,blob in extra.items():manifest[name]=hashlib.sha256(blob).hexdigest();z.writestr(name,blob)
    z.writestr('SHA256SUMS.json',json.dumps(manifest,indent=2))
with zipfile.ZipFile(temp) as z:
    assert z.testzip() is None
    for name,sha in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==sha
    assert hashlib.sha256(z.read('inputs/original_collaborator.cif')).hexdigest()==summary['source_cif_sha256']
temp.replace(archive)
shutil.copy2(out/'report.html',results/'collaborator_CIF_review.html')
shutil.copy2(out/'ed_scores.csv',results/'collaborator_ED_scores.csv')
records=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
         for p in sorted(results.iterdir()) if p.is_file() and p.name!='manifest.json']
(results/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'hashed_files':len(manifest)},indent=2))
