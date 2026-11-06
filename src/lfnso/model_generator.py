"""Hypothesis structures from supplied coordinates and user-defined Fe-site Ni substitution."""
import itertools
import numpy as np
from scipy.optimize import linear_sum_assignment
from pymatgen.core import Structure,Lattice,Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

FE_NI={'Fe':.67,'Ni':.33}


def sb_orbits(template):
    ids=[i for i,s in enumerate(template) if s.species_string=='Mn']
    if len(ids)!=4:raise ValueError('Expected four Mn sites in the supplied conventional Li2MnO3 cell')
    reps=set()
    for combo in itertools.combinations(ids,2):
        orbit=[]
        for op in SpacegroupAnalyzer(template,symprec=.01).get_symmetry_operations():
            mapped=[]
            for i in combo:
                dist=template.lattice.get_all_distances([op.operate(template[i].frac_coords)],[template[j].frac_coords for j in ids])[0]
                if min(dist)>.01:raise ValueError('Symmetry mapping failed')
                mapped.append(ids[int(np.argmin(dist))])
            orbit.append(tuple(sorted(mapped)))
        reps.add(min(orbit))
    return sorted(reps)


def make_parent(template,sb_indices,lattice_parameters):
    species=[]
    for i,site in enumerate(template):
        species.append('Sb' if i in sb_indices else FE_NI if site.species_string=='Mn' else site.species)
    return Structure(Lattice.from_parameters(*lattice_parameters),species,template.frac_coords,
                     site_properties={'parent_index':list(range(len(template))),
                                      'parent_y':template.frac_coords[:,1].tolist(),
                                      'parent_z':template.frac_coords[:,2].tolist()})


def generate(parent,ny=1,nz=1,harmonic=1,exchange=0.,displacement_A=0.,stacking=False,seed=0,integer_approx=False,host_exchange=0.):
    if ny<1 or nz<1 or int(ny)!=ny or int(nz)!=nz:raise ValueError('Positive integer dimensions required')
    if not 0<=exchange<=1:raise ValueError('Exchange amplitude must be between zero and one')
    s=parent.copy();s.make_supercell([1,ny,nz])
    # Preserve original slot-to-tile mapping before any displacement.
    slots={}
    for i,site in enumerate(s):
        j=int(round(site.frac_coords[1]*ny-site.properties['parent_y']))%ny
        k=int(round(site.frac_coords[2]*nz-site.properties['parent_z']))%nz
        slots[(site.properties['parent_index'],j,k)]=i
    if exchange:
        li=[i for i,t in enumerate(parent) if t.species_string=='Li' and abs(t.frac_coords[2])<1e-6]
        fe=[i for i,t in enumerate(parent) if 'Fe' in t.species]
        if len(li)!=len(fe):raise ValueError('Need matched in-slab Li and Fe-host sites')
        distance=parent.lattice.get_all_distances(parent.frac_coords[li],parent.frac_coords[fe])
        aa,bb=linear_sum_assignment(distance)
        for j in range(ny):
            amount=exchange*(1+np.cos(2*np.pi*harmonic*j/ny))/2
            for k in range(nz):
                for a,b in zip(aa,bb):
                    i=slots[(li[a],j,k)]
                    s.replace(i,{'Li':1-amount,'Fe':.67*amount,'Ni':.33*amount},properties=s[i].properties.copy())
                    i=slots[(fe[b],j,k)]
                    s.replace(i,{'Li':amount,'Fe':.67*(1-amount),'Ni':.33*(1-amount)},properties=s[i].properties.copy())
    if host_exchange:
        if exchange or not 0<=host_exchange<=1:raise ValueError('Use a separate bounded Sb/Fe-host exchange control')
        fe=[i for i,t in enumerate(parent) if 'Fe' in t.species]
        sb=[i for i,t in enumerate(parent) if t.species_string=='Sb']
        for j in range(ny):
            amount=host_exchange*(1+np.cos(2*np.pi*harmonic*j/ny))/2
            for k in range(nz):
                for a,b in zip(fe,sb):
                    i=slots[(a,j,k)];s.replace(i,{'Sb':amount,'Fe':.67*(1-amount),'Ni':.33*(1-amount)},properties=s[i].properties.copy())
                    i=slots[(b,j,k)];s.replace(i,{'Sb':1-amount,'Fe':.67*amount,'Ni':.33*amount},properties=s[i].properties.copy())
    if integer_approx:
        if exchange:raise ValueError('Integer approximation implemented only on unexchanged Fe hosts')
        ids=[i for i,t in enumerate(s) if 'Fe' in t.species]
        if len(ids)%3:raise ValueError('2/3 Fe approximant requires Fe-host count divisible by 3')
        chosen=np.random.default_rng(seed).permutation(ids);ni=set(chosen[:len(ids)//3])
        for i in ids:s.replace(i,'Ni' if i in ni else 'Fe',properties=s[i].properties.copy())
    if displacement_A:
        # A coherent all-atom acoustic trial; not an energy-relaxed oxygen structure.
        direction=parent.lattice.matrix[1]/parent.lattice.b
        for i,t in enumerate(s):
            u=displacement_A*np.cos(2*np.pi*harmonic*t.frac_coords[1])
            s.translate_sites([i],u*direction,frac_coords=False,to_unit_cell=True)
    if stacking:
        if nz!=4:raise ValueError('This AABB trial requires four slabs')
        shifts=np.array([0,0,1/6,1/6])
        for i,t in enumerate(s):
            # Assign O on either side to its closest TM slab; translate slab Li/TM/O together.
            z=t.frac_coords[2]*nz
            interlayer=t.properties['parent_z']==.5 and 'Li' in t.species
            if interlayer:
                lo=int(np.floor(z))%nz;amount=(shifts[lo]+shifts[(lo+1)%nz])/2
            else:amount=shifts[int(np.floor(z+.5))%nz]
            s.translate_sites([i],amount*parent.lattice.matrix[1],frac_coords=False,to_unit_cell=True)
    return s


def validate(s,formula_units):
    expected={'Li':4.,'Fe':.67,'Ni':.33,'Sb':1.,'O':6.}
    normalized={k:float(v/formula_units) for k,v in s.composition.as_dict().items()}
    exact=all(abs(normalized.get(k,0)-v)<1e-7 for k,v in expected.items())
    oxygen=[i for i,t in enumerate(s) if t.species_string=='O']
    cations=[i for i in range(len(s)) if i not in oxygen]
    distances=s.distance_matrix.copy();np.fill_diagonal(distances,np.inf)
    coordination=[]
    for i in cations:
        ds=sorted(n.nn_distance for n in s.get_neighbors(s[i],2.7) if n.species_string=='O')
        coordination.append({'site':i,'species':s[i].species_string,'oxygen_neighbors_below_2p7_A':len(ds),'distances_A':list(map(float,ds))})
    return {'composition_per_Li4_formula':normalized,'exact_decimal_composition':exact,
            'all_site_occupancy_sums_one':all(abs(sum(t.species.values())-1)<1e-9 for t in s),
            'minimum_distinct_site_distance_A':float(distances.min()),
            'sixfold_O_coordination_fraction':sum(t['oxygen_neighbors_below_2p7_A']==6 for t in coordination)/len(coordination),
            'coordination':coordination,'coordination_cutoff_A':2.7,
            'gross_geometry_pass':bool(distances.min()>1.4),
            'chemical_status':'unrelaxed inherited coordinates; coordination/distance checks are necessary, not energetic validation',
            'charge_per_formula_scenarios':{'Li1_Fe3_Ni2_Sb5_Ominus2':normalized.get('Li',0)+3*normalized.get('Fe',0)+2*normalized.get('Ni',0)+5*normalized.get('Sb',0)-2*normalized.get('O',0),'Li1_Fe3_Ni3_Sb5_Ominus2':normalized.get('Li',0)+3*normalized.get('Fe',0)+3*normalized.get('Ni',0)+5*normalized.get('Sb',0)-2*normalized.get('O',0)},
            'charge_note':'Scenarios only; oxidation states and oxygen nonstoichiometry are not inferred or imposed.'}
