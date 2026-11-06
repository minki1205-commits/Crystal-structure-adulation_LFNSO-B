"""Separate c-layer occupation, registry, breathing and coupled hypotheses."""
import argparse,base64,csv,html,json,time
from pathlib import Path
import numpy as np
from pymatgen.core import Lattice
from pymatgen.io.cif import CifParser
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from scipy.spatial import cKDTree
from .coupled_models import extend_layers,primitive_binary_necklaces,plane_reflections_3d
from .occupational_search import form_factor,ProfileObjective
from .model_generator import validate
from .ed_matching import electron_intensity,transform,in_view,match_one_to_one
from .ed_metric_diagnostic import apply_correction
from .pipeline import digest,write_json


def powder(structure,wavelength=1.5406):
    lo,hi=2*np.sin(np.deg2rad(np.array([16,46])/2))/wavelength
    pts=structure.lattice.reciprocal_lattice_crystallographic.get_points_in_sphere([[0,0,0]],[0,0,0],hi)
    pairs=[(np.rint(p[0]).astype(int),p[1]) for p in pts if p[1]>=lo and tuple(np.rint(p[0]).astype(int))>(0,0,0)]
    pairs.sort(key=lambda p:(p[1],*-p[0]));hk=np.array([p[0] for p in pairs]);g=np.array([p[1] for p in pairs])
    theta=np.arcsin(wavelength*g/2);angle=np.rad2deg(theta)*2
    _,first,groups=np.unique(np.rint(angle/XRDCalculator.TWO_THETA_TOL).astype(np.int64),return_index=True,return_inverse=True)
    fs={e.symbol:form_factor(e.symbol,g) for e in structure.elements};intensity=[]
    for ids in np.array_split(np.arange(len(hk)),max(1,int(np.ceil(len(hk)/256)))):
        phase=np.exp(2j*np.pi*(hk[ids]@structure.frac_coords.T));f=np.zeros(len(ids),complex)
        for i,site in enumerate(structure):
            for element,amount in site.species.items():f+=phase[:,i]*fs[element.symbol][ids]*amount
        intensity.extend(abs(f)**2)
    lp=2*(1+np.cos(2*theta)**2)/(np.sin(theta)**2*np.cos(theta))
    return angle[first],np.bincount(groups,weights=np.array(intensity)*lp,minlength=len(first))


def verify_powder(structure,angles,intensities,wavelength):
    expected=XRDCalculator(wavelength=wavelength).get_pattern(structure,two_theta_range=(16,46))
    relative=intensities/intensities.max()*100;mask=relative>XRDCalculator.SCALED_INTENSITY_TOL
    np.testing.assert_allclose(angles[mask],expected.x,atol=1e-7)
    np.testing.assert_allclose(relative[mask],expected.y,atol=2e-6,rtol=2e-6)
    return {'peak_count':len(expected.x),'maximum_relative_intensity_difference':float(np.max(abs(relative[mask]-expected.y)))}


def specifications():
    specs=[{'id':'base_1c','family':'base','sequence':[0],'operation':'none'},
           {'id':'exact_repeat_2c','family':'repeat_control','sequence':[0,0],'operation':'none'}]
    for n in [2,3,4,6]:
        for sequence in primitive_binary_necklaces(n):
            specs.append({'id':'SbFe_'+''.join(map(str,sequence))+f'_{n}c','family':'c_occupation','sequence':list(sequence),'operation':'swap_Sb_Fe'})
    specs.extend([
        {'id':'LiSbFe_ABC_3c','family':'c_occupation_Li','sequence':[0,1,2],'operation':'cycle_Li_Sb_Fe'},
        {'id':'LiSbFe_AABBCC_6c','family':'c_occupation_Li','sequence':[0,0,1,1,2,2],'operation':'cycle_Li_Sb_Fe'},
        *[{'id':f'registry_AABB_b{denominator}_4c','family':'rigid_registry','sequence':[0,0,1,1],'operation':'none','registry_fraction':1/denominator} for denominator in [48,24,6]],
        {'id':'breathing_3c_004','family':'cation_breathing','sequence':[0,0,0],'operation':'none','breathing_A':.04},
        {'id':'breathing_6c_008','family':'cation_breathing','sequence':[0]*6,'operation':'none','breathing_A':.08},
        {'id':'coupled_SbFe_AABB_b48_4c','family':'coupled','sequence':[0,0,1,1],'operation':'swap_Sb_Fe','registry_fraction':1/48},
    ])
    return specs


def run(config_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cp=Path(config_path).resolve();cfg=json.loads(cp.read_text());out=(cp.parent/cfg['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    ed_path=(cp.parent/cfg['ed_summary']).resolve();affine_path=(cp.parent/cfg['affine_summary']).resolve()
    ed=json.loads(ed_path.read_text());affine=json.loads(affine_path.read_text());data=np.loadtxt(cfg['xrd_csv'],delimiter=',',skiprows=1)
    write_json(out/'provenance.json',{'config':cfg,'code_sha256':{name:digest(Path(__file__).with_name(name)) for name in ['coupled_models.py','coupled_search.py','ed_matching.py','ed_metric_diagnostic.py','occupational_search.py','model_generator.py']},'ed_summary_sha256':digest(ed_path),'affine_summary_sha256':digest(affine_path),'xrd_sha256':digest(Path(cfg['xrd_csv']))})
    # Geometry-only diagnostic: more cells add allowed positions, not evidence for a period.
    grid_rows=[]
    for dataset in ed['datasets']:
        if dataset['nominal_zone']!=[1,1,0]:continue
        key=dataset['id']+'__conditional_B__1,1,0';reg=ed['registrations'][key];correction=affine['corrections'][key]
        lattice=Lattice.from_parameters(*cfg['bases'][-1]['parent_metric'])
        for n in [1,2,3,4,6,7,8,12]:
            hk,hp,g,xy=plane_reflections_3d(lattice,[1,1,n],[1,1,0],14)
            pred=apply_correction(transform(xy,reg['parameters'],reg['mirror']),correction);pred=pred[in_view(pred,dataset)]
            match=match_one_to_one(dataset['observed'],pred,dataset['position_tolerances_nm_inv'][1]);near=cKDTree(pred).query(dataset['observed'])[0]
            grid_rows.append({'dataset':dataset['id'],'c_repeat':n,'median_nearest_distance_nm_inv':float(np.median(near)),**{k:v for k,v in match.items() if k!='pairs'}})
    write_json(out/'geometry_only_c_grid.json',grid_rows)
    rows=[];t0=time.monotonic();windows={'superstructure':[20.4,34.6],'high_angle':[35,45],'low_angle':[17,20.3]}
    for base_entry in cfg['bases']:
        base_path=(cp.parent/base_entry['cif']).resolve();base=CifParser(base_path,frac_tolerance=0).parse_structures(primitive=False)[0]
        ny=base_entry['b_repeat'];lattice=Lattice.from_parameters(*base_entry['parent_metric']);base_profile=None
        for spec in specifications():
            identifier=base_entry['id']+'__'+spec['id'];folder=out/identifier;folder.mkdir()
            args={k:v for k,v in spec.items() if k not in ['id','family']};s=extend_layers(base,base_b_repeat=ny,**args);nz=len(spec['sequence']);checks=validate(s,2*ny*nz)
            scores={};profiles={};angles,intensity=powder(s,cfg['wavelength_A'])
            for name,(lo,hi) in windows.items():
                mask=(data[:,0]>=lo)&(data[:,0]<=hi);obj=ProfileObjective(angles,*data[mask].T,cfg['width_deg']);v=obj.evaluate(intensity)
                profiles[name]=v['profile'];scores[name]={k:value for k,value in v.items() if k!='profile'}
            control=None
            if spec['family']=='base':base_profile=np.array(profiles['superstructure']);control=verify_powder(s,angles,intensity,cfg['wavelength_A'])
            if spec['family']=='repeat_control':
                np.testing.assert_allclose(profiles['superstructure'],base_profile,atol=1e-6,rtol=1e-8);control={'exact_repeat_profile_matches_base':True}
            if spec['id']=='SbFe_000001_6c' and base_entry['id']==cfg['bases'][-1]['id']:control=verify_powder(s,angles,intensity,cfg['wavelength_A'])
            s.to(filename=str(folder/'model.cif'),fmt='cif',symprec=None)
            write_json(folder/'xrd.json',{'two_theta_deg':angles.tolist(),'intensity_arbitrary_common_units':intensity.tolist(),'wavelength_A':cfg['wavelength_A']})
            edscores=[]
            for dataset in ed['datasets']:
                zone=[1,0,-1] if dataset['id']=='DP01_101' else dataset['nominal_zone']
                key=dataset['id']+'__'+base_entry['metric_name']+'__'+','.join(map(str,zone));reg=ed['registrations'][key]
                hk,hp,g,xy=plane_reflections_3d(lattice,[1,ny,nz],zone,14)
                values=electron_intensity(s,hk,np.linalg.norm(g,axis=1));prediction=transform(xy,reg['parameters'],reg['mirror'])
                # Store reflection tables so assignments remain inspectable; image pose is never re-fitted.
                write_json(folder/(dataset['id']+'_reflections.json'),{'zone_hypothesis':zone,'hkl_supercell':hk.tolist(),'hkl_parent':hp.tolist(),'g_similarity_nm_inv':prediction.tolist(),'kinematic_intensity':values.tolist()})
                for mode in ['similarity','affine_diagnostic']:
                    pp=prediction if mode=='similarity' else apply_correction(prediction,affine['corrections'][key]);vis=in_view(pp,dataset)
                    norm=max(float(values[vis].max()) if vis.any() else 0,1e-30);relative=100*values/norm
                    for threshold in [.1,1,5]:
                        selected=vis&(relative>=threshold);ids=np.flatnonzero(selected)
                        for tolerance in dataset['position_tolerances_nm_inv']:
                            m=match_one_to_one(dataset['observed'],pp[selected],tolerance)
                            for pair in m['pairs']:pair['pattern_reflection_index']=int(ids[pair['predicted_index']]);pair['source_peak_index']=dataset['source_peak_indices'][pair['observed_index']]
                            edscores.append({'dataset':dataset['id'],'mode':mode,'zone_hypothesis':zone,'threshold_percent':threshold,'tolerance_nm_inv':tolerance,**m})
            metadata={'id':identifier,'base_cif_sha256':digest(base_path),'base':base_entry,'specification':spec,'site_count':len(s),'repeats':[1,ny,nz],'validation':checks,'xrd_scores':scores,'control_validation':control,'ed_scores':edscores,
                      'geometry_eligible':bool(checks['exact_decimal_composition'] and checks['all_site_occupancy_sums_one'] and checks['gross_geometry_pass'] and checks['sixfold_O_coordination_fraction']==1),
                      'interpretation':'Unrelaxed hypothesis. Sequence enumeration identifies cyclic sequence classes, not proven space-group-distinct 3D structures. ED pose reused from parent; no candidate-specific calibration.',
                      'gate1':'not established','gate2':'not established','STEM':'not run'}
            write_json(folder/'metadata.json',metadata)
            row={k:v for k,v in metadata.items() if k not in ['validation','ed_scores']};row['geometry_summary']={k:v for k,v in checks.items() if k!='coordination'}
            row['ed_standard']=[{k:v for k,v in m.items() if k!='pairs'} for m in edscores if m['threshold_percent']==1 and m['tolerance_nm_inv']==(.08 if m['dataset'] in ['DP01_101','dp3_110'] else .24)]
            rows.append(row)
            print(identifier,'XRD',round(scores['superstructure']['R2_vs_linear_background'],4),'CN6',round(checks['sixfold_O_coordination_fraction'],4),'110 F1',[(r['mode'],round(r['f1'],3)) for r in row['ed_standard'] if r['dataset']=='dp3_110'],flush=True)
    summary={'config':cfg,'elapsed_seconds':time.monotonic()-t0,'candidate_count':len(rows),'geometry_eligible_count':sum(r['geometry_eligible'] for r in rows),'rows':rows,'validated_shortlist':[],
             'limitations':['This adaptive exploratory extension uses the previously inspected XRD and ED data; no fresh validation set.',
                            'Cation species are exact decimal average occupancies. Oxidation states, oxygen vacancies and energy stability are not assigned.',
                            'Powder profile uses provisional Cu Kalpha and one fixed Gaussian width; no full refinement.',
                            'Binary layer sequences quotient cyclic translation only; reversal or chemical swaps are not assumed equivalent.',
                            'c-grid coverage alone grows with grid density and does not select a period.',
                            'Registry moves slab O together with cations and averages neighboring shifts for interlayer Li. Large shifts can violate coordination and are excluded from eligible comparison.',
                            'Breathing amplitude is the Fe/Sb displacement amplitude; slab Li moves at -2 times it. Oxygen is fixed for this separate displacement test.',
                            'No dynamical ED, tilt, ionic scattering or STEM simulation; existing axis/calibration ambiguities persist.']}
    write_json(out/'summary.json',summary)
    with (out/'scores.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['id','family','geometry_eligible','XRD_R2_superstructure','XRD_R2_high_angle','dataset','mode','ED_F1','ED_matches','ED_predicted'])
        for r in rows:
            for m in r['ed_standard']:w.writerow([r['id'],r['specification']['family'],r['geometry_eligible'],r['xrd_scores']['superstructure']['R2_vs_linear_background'],r['xrd_scores']['high_angle']['R2_vs_linear_background'],m['dataset'],m['mode'],m['f1'],m['matches'],m['predicted_count']])
    page='<!doctype html><meta charset="utf-8"><title>LFNSO b/c coupled hypotheses</title><style>body{font:16px system-ui;max-width:1300px;margin:30px auto}td,th{padding:6px;border-bottom:1px solid #ccc}img{max-width:100%}pre{white-space:pre-wrap}</style><h1>c축 점유·적층·변위 가설 확장</h1><p>기존 7b/11b 후보 세 개에서 출발해 c축의 점유 sequence, 작은 registry shift, slab cation breathing과 결합 가설을 구분했습니다. 조성·배위 검사와 기존 pose의 다중 ED 위치 비교를 수행했습니다. 최종 구조가 아니며 gate 통과를 선언하지 않습니다.</p>'
    page+='<h2>기하 검사 통과 후보 · XRD 점수 순서, 공동 순위 아님</h2><table><tr><th>후보</th><th>XRD super R²</th><th>high R²</th><th>[110] DM3 F1 / affine diagnostic</th></tr>'
    for r in sorted([r for r in rows if r['geometry_eligible']],key=lambda r:r['xrd_scores']['superstructure']['R2_vs_linear_background'],reverse=True):
        edv=next(m['f1'] for m in r['ed_standard'] if m['dataset']=='dp3_110' and m['mode']=='affine_diagnostic')
        page+=f'<tr><td>{r["id"]}</td><td>{r["xrd_scores"]["superstructure"]["R2_vs_linear_background"]:.4f}</td><td>{r["xrd_scores"]["high_angle"]["R2_vs_linear_background"]:.4f}</td><td>{edv:.4f}</td></tr>'
    page+='</table><pre>'+html.escape(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre>'
    (out/'report.html').write_text(page)
    snapshots=out/'source_snapshot';snapshots.mkdir()
    for name in ['coupled_models.py','coupled_search.py']:(snapshots/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(p.parse_args().config)
