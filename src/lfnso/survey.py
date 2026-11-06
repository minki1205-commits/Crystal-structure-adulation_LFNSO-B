"""Run an auditable first survey on actual CIF/DM inputs."""
import argparse,csv,html,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from .pipeline import write_json,digest
from .data_io import read_dm
from .experimental import calibrated_fft,extract_peaks


def preview(a,path,title,peaks=None,reciprocal=False):
    with np.errstate(invalid='ignore'):
        a=np.abs(a).astype(float)
    finite=np.isfinite(a)
    if reciprocal:a=np.log1p(np.maximum(a,0))
    lo,hi=np.percentile(a[finite],[1,99.8])
    shown=np.ma.masked_where(~finite,a)
    fig,ax=plt.subplots(figsize=(7,7))
    ax.imshow(shown,cmap='gray',vmin=lo,vmax=hi)
    if peaks:
        ps=peaks['peaks'];ax.scatter([p['x_px'] for p in ps],[p['y_px'] for p in ps],s=20,facecolors='none',edgecolors='tomato',linewidths=.5)
        c=peaks['center_yx_px'];ax.plot(c[1],c[0],'+',color='cyan')
    ax.set(title=title,xlabel='x pixel',ylabel='y pixel');fig.tight_layout();fig.savefig(path,dpi=130);plt.close(fig)


def run(config_path):
    config_path=Path(config_path).resolve();cfg=json.loads(config_path.read_text())
    root=(config_path.parent/cfg['input_directory']).resolve();out=(config_path.parent/cfg['output']).resolve()
    if out==root or out.is_relative_to(root) or root.is_relative_to(out):raise ValueError('Separate input/output trees required')
    out.mkdir(parents=True,exist_ok=False)
    for n in ['previews','peaks','ffts']: (out/n).mkdir()
    write_json(out/'config.json',cfg)
    records=[];arrays={};sections=[]
    for p in sorted(root.glob('*.dm*')):
        a,row=read_dm(p);arrays[p.name]=a
        zone=re.search(r'_(\d{3})\.dm[34]$',p.name)
        row.update({'filename':p.name,'zone_axis_from_filename':zone.group(1) if zone else None,
                    'zone_axis_validation':'not independently indexed',
                    'sample_identity':cfg.get('sample_identity','unconfirmed'),
                    'modality':'stored_image_fft' if p.name.startswith('FFT of ') else ('diffraction' if row['pixel_unit_yx']==['1/nm','1/nm'] else 'real_space_image'),
                    'detector_identity':('BF' if p.name.startswith('BF_') else cfg.get('image1_detector','unconfirmed') if 'image1' in p.name else 'diffraction camera'),
                    'identity_provenance':cfg.get('identity_provenance','unconfirmed'),
                    'acquisition_group':p.stem.split('_')[-2] if 'image1' in p.name or p.name.startswith('BF_') else p.stem})
        key=re.sub(r'[^A-Za-z0-9_-]','_',p.stem)
        peaks=None
        if row['quantitative_read_status']=='readable':
            if row['pixel_unit_yx']==['1/nm','1/nm']:
                peaks=extract_peaks(a,row['pixel_size_yx'],is_fft=row['modality']=='stored_image_fft')
            elif row['pixel_unit_yx']==['nm','nm']:
                row['field_of_view_yx_nm']=(np.array(a.shape)*row['pixel_size_yx']).tolist()
                fft,delta=calibrated_fft(a,row['pixel_size_yx'])
                row['fft_bin_yx_nm_inv']=delta.tolist()
                row['fft_processing']='mean subtraction + separable Hann window; nominal bin size is not peak uncertainty'
                fp=extract_peaks(fft,delta,is_fft=True)
                write_json(out/'ffts'/f'{key}.json',fp)
                preview(fft,out/'previews'/f'{key}_fft.png',p.name+' / Hann FFT',fp,True)
                row['fft_preview']=f'previews/{key}_fft.png'
        if peaks:
            write_json(out/'peaks'/f'{key}.json',peaks)
            with (out/'peaks'/f'{key}.csv').open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['y_px','x_px','g_y_nm-1','g_x_nm-1','g_norm_nm-1','classification'])
                for q in peaks['peaks']:w.writerow([q['y_px'],q['x_px'],*q['g_yx_nm_inv'],q['g_norm_nm_inv'],q['classification']])
            row['peak_count']=len(peaks['peaks']);row['center_summary']={k:peaks[k] for k in ['center_yx_px','center_method','opposite_pairs','midpoint_rms_yx_px']}
        preview(a,out/'previews'/f'{key}.png',p.name,peaks,row['modality']!='real_space_image')
        row['preview']=f'previews/{key}.png';records.append(row)
    original='image1_20MX_200kV_0075_110.dm4';stored='FFT of image1_20MX_200kV_0075_110.dm3'
    fft_check=None
    if original in arrays and stored in arrays:
        calculated=np.fft.fftshift(np.fft.fft2(arrays[original]));saved=arrays[stored]
        mask=np.ones(saved.shape,bool);mask[saved.shape[0]//2,saved.shape[1]//2]=False
        fft_check={'source':original,'derived':stored,'comparison':'unwindowed full-frame FFT magnitude, DC excluded',
                   'relative_magnitude_l2_error_excluding_DC':float(np.linalg.norm((np.abs(saved)-np.abs(calculated))[mask])/np.linalg.norm(np.abs(calculated)[mask])),
                   'complex_phase_convention':'not established; complex arrays differ',
                   'independent_constraint':False}
    cif=root/cfg['reference_cif'];parser=CifParser(cif);s=parser.parse_structures(primitive=False)[0]
    symmetry=SpacegroupAnalyzer(s,symprec=.01);sym=symmetry.get_symmetrized_structure()
    cif_report={'source':str(cif),'sha256':digest(cif),'role':'Li2MnO3 geometric reference, NOT LFNSO chemical parent',
                'composition':s.composition.as_dict(),'reduced_formula':s.composition.reduced_formula,'site_count':len(s),
                'cell_parameters':list(s.lattice.parameters),'volume_A3':s.volume,
                'detected_space_group':symmetry.get_space_group_symbol(),'symprec_A':.01,
                'site_groups':[{'wyckoff':w,'species':str(s[g[0]].species),'multiplicity':len(g),'representative_frac':s[g[0]].frac_coords.tolist()} for w,g in zip(sym.wyckoff_symbols,sym.equivalent_indices)],
                'occupancy_sum_range':[float(min(sum(site.species.values()) for site in s)),float(max(sum(site.species.values()) for site in s))],
                'nominal_charge_check':{'assumed_valences':{'Li':1,'Mn':4,'O':-2},'cell_charge_e':float(s.composition['Li']+4*s.composition['Mn']-2*s.composition['O']),'status':'conditional on conventional valences, not measured'},
                'lfnso_prompt_metric_only':[5.1706,8.9382,5.1635,90,109.49,90]}
    q=[]
    for period in [3,4,6,12]:
        q.append({'period_b':period,'reference_q_nm_inv':10/(s.lattice.b*period),'lfnso_prompt_q_nm_inv':10/(8.9382*period),'status':'hypothesis, not measured'})
    result={'status':'input_survey_complete_structural_ranking_deferred','files':records,'cif':cif_report,'fft_dependency_check':fft_check,
            'q_hypotheses':q,'unavailable':cfg.get('unavailable',[]),'xrd':'not supplied','ranked_candidates':[],
            'gate1':'not_evaluated','gate2':'not_evaluated','stem_simulation':'not_run',
            'missing':['LFNSO parent CIF/site assignment','experimental XRD and wavelength','sample identity per acquisition','HAADF detector inner/outer angles, convergence angle and thickness range','independently validated hkl/orientation/calibration']}
    write_json(out/'survey.json',result)
    table=''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['filename'],r['shape_yx'],r['pixel_size_yx'],r['pixel_unit_yx'],r['quantitative_read_status'],r.get('peak_count','—')])+'</tr>' for r in records)
    for r in records:
        sections.append('<h2>'+html.escape(r['filename'])+'</h2><img src="'+r['preview']+'">'+('<img src="'+r['fft_preview']+'">' if 'fft_preview' in r else ''))
    page='''<!doctype html><html lang="ko"><meta charset="utf-8"><title>LFNSO experimental input survey</title>
<style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:20px}td,th{border:1px solid #ccc;padding:8px}table{border-collapse:collapse}img{max-width:48%}pre{white-space:pre-wrap}</style>
<h1>LFNSO21-B 입력 자료 1차 점검</h1><p>시료 LFNSO21-B 및 image1=HAADF는 사용자 확인 정보입니다. Zone axis는 파일명 기준으로, 아직 독립 지수화되지 않았습니다. 점 표시는 자동 검출 후보이며 hkl/satellite 확정이 아닙니다. 원본은 보존했고 비정상 값을 가진 OneView는 정량 비교에서 제외했습니다.</p>
<p>Li₂MnO₃ CIF는 기하학적 기준입니다. LFNSO 조성으로 치환하지 않았고 XRD/ED gate 및 구조 순위를 산출하지 않았습니다.</p>'''
    page+='<table><tr><th>파일</th><th>크기 y,x</th><th>pixel size y,x</th><th>단위</th><th>상태</th><th>검출 후보</th></tr>'+table+'</table>'
    page+='<h2>CIF / 데이터 중복 검증 / 가설 q</h2><pre>'+html.escape(json.dumps({'cif':cif_report,'fft_check':fft_check,'q':q,'unavailable':result['unavailable'],'missing':result['missing']},indent=2,ensure_ascii=False))+'</pre>'+''.join(sections)+'</html>'
    (out/'report.html').write_text(page)
    print(json.dumps({'output':str(out),'files':len(records),'status':result['status']}))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
