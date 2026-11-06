"""Explicit b/c supercells and independent occupation, registry, displacement trials."""
import itertools
import numpy as np
from pymatgen.core import Structure
from .model_generator import FE_NI


def plane_reflections_3d(parent_lattice, repeats, zone, gmax):
    dims=np.asarray(repeats,int);z=np.asarray(zone,int)
    if np.any(dims<1) or not np.any(z):raise ValueError('Positive repeats and nonzero zone required')
    common=np.lcm.reduce(dims);coefficient=z*(common//dims)
    bounds=np.ceil(np.array(parent_lattice.abc)*dims*gmax/10).astype(int)+1
    solved=int(np.flatnonzero(coefficient)[-1]);free=[i for i in range(3) if i!=solved]
    aa,bb=np.meshgrid(np.arange(-bounds[free[0]],bounds[free[0]]+1),np.arange(-bounds[free[1]],bounds[free[1]]+1),indexing='ij')
    hk=np.zeros((aa.size,3),int);hk[:,free[0]]=aa.ravel();hk[:,free[1]]=bb.ravel();num=-(hk@coefficient)
    mask=num%coefficient[solved]==0;hk=hk[mask];hk[:,solved]=num[mask]//coefficient[solved]
    hp=hk/dims;g=hp@parent_lattice.reciprocal_lattice_crystallographic.matrix*10
    norm=np.linalg.norm(g,axis=1);mask=(norm>1e-8)&(norm<=gmax+1e-9);hk,hp,g=hk[mask],hp[mask],g[mask]
    normal=z@parent_lattice.matrix;normal/=np.linalg.norm(normal)
    seed=np.eye(3)[np.argmin(abs(normal))];e1=np.cross(normal,seed);e1/=np.linalg.norm(e1);e2=np.cross(normal,e1)
    return hk,hp,g,g@np.array([e1,e2]).T


def primitive_binary_necklaces(period):
    """Cyclic translations only. Do not quotient reversal or swap chemical identities."""
    result=[]
    for seq in itertools.product([0,1],repeat=period):
        if len(set(seq))<2:continue
        if any(period%d==0 and seq==seq[:d]*(period//d) for d in range(1,period)):continue
        if seq!=min(seq[i:]+seq[:i] for i in range(period)):continue
        result.append(seq)
    return result


def extend_layers(base, sequence, operation='swap_Sb_Fe', registry_fraction=0., breathing_A=0., base_b_repeat=1):
    """Preserve counts; Ni remains .67/.33 with Fe. No energy relaxation implied."""
    sequence=np.asarray(sequence,int);nz=len(sequence)
    if not nz:raise ValueError('Layer sequence required')
    p=base.copy();p.add_site_property('base_z',np.mod(p.frac_coords[:,2],1).tolist())
    s=p.copy();s.make_supercell([1,1,nz])
    normal=np.cross(base.lattice.matrix[0],base.lattice.matrix[1]);normal/=np.linalg.norm(normal)
    shifts=sequence.astype(float)*registry_fraction
    for i,site in enumerate(s):
        z=site.frac_coords[2]*nz
        layer=int(np.floor(z+.5))%nz
        slab=abs(site.properties['base_z']-np.rint(site.properties['base_z']))<1e-7 and site.species_string!='O'
        state=int(sequence[layer])
        if slab and operation=='swap_Sb_Fe' and state%2:
            if 'Sb' in site.species:s.replace(i,FE_NI,properties=site.properties.copy())
            elif 'Fe' in site.species:s.replace(i,'Sb',properties=site.properties.copy())
        elif slab and operation=='cycle_Li_Sb_Fe':
            label=0 if site.species_string=='Li' else 1 if 'Sb' in site.species else 2
            s.replace(i,['Li','Sb',FE_NI][(label+state)%3],properties=site.properties.copy())
        elif operation not in ['none','swap_Sb_Fe','cycle_Li_Sb_Fe']:
            raise ValueError('Unknown layer operation')
        if registry_fraction:
            interlayer=abs(site.properties['base_z']-.5)<1e-7 and site.species_string=='Li'
            if interlayer:
                lower=int(np.floor(z))%nz;fraction=(shifts[lower]+shifts[(lower+1)%nz])/2
            else:fraction=shifts[layer]
            s.translate_sites([i],fraction*base.lattice.matrix[1]/base_b_repeat,frac_coords=False,to_unit_cell=True)
        if breathing_A and slab:
            weight=-2. if s[i].species_string=='Li' else 1.
            displacement=breathing_A*np.cos(2*np.pi*layer/nz)*weight
            s.translate_sites([i],displacement*normal,frac_coords=False,to_unit_cell=True)
    for site in s:site.label=site.species_string
    return s
