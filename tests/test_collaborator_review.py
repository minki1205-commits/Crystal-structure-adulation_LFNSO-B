"""Independent metric identities and disordered-CIF export regression checks."""
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifParser
from lfnso.collaborator_review import predict_positions, metric_fit, export_calculation_copy


def test_monoclinic_positions_against_analytic_reciprocal_metric():
    a,b,c,beta=5.16,8.94,9.87,100.03
    hk=np.array([[0,0,2],[1,1,0],[0,2,1],[1,1,-2],[1,1,2],[0,2,3]])
    h,k,l=hk.T;angle=np.deg2rad(beta)
    inverse_d2=(h*h/a**2+l*l/c**2-2*h*l*np.cos(angle)/(a*c))/np.sin(angle)**2+k*k/b**2
    analytic=np.rad2deg(2*np.arcsin(1.5406*np.sqrt(inverse_d2)/2))
    np.testing.assert_allclose(predict_positions([a,b,c,beta],hk,1.5406),analytic,atol=1e-12)


def test_metric_fit_recovers_known_cell_preserving_coordinates_and_occupancy():
    structure=Structure(Lattice.monoclinic(5.16,8.94,9.87,100.03),[{'Li':.5,'Fe':1/6,'Ni':1/3},'O'],[[.1,.2,.3],[.2,.4,.6]])
    target=[5.28,8.78,9.91,100.77]
    hk=np.array([[0,0,2],[1,1,0],[0,2,1],[1,1,-2],[1,1,2],[0,2,3]])
    observed=predict_positions(target,hk,1.5406)
    fit,report=metric_fit(structure,observed,hk,1.5406)
    np.testing.assert_allclose(report['parameters_a_b_c_beta'],target,atol=1e-6)
    np.testing.assert_allclose(fit.frac_coords,structure.frac_coords)
    assert fit.composition==structure.composition


def test_export_has_unique_element_labels_and_preserves_disordered_sites(tmp_path):
    structure=Structure(Lattice.monoclinic(5.16,8.94,9.87,100.03),[{'Li':.5,'Fe':.1667,'Ni':.3333},'O'],[[.1,.2,.3],[.2,.4,.6]],labels=['mixed','oxygen'])
    path=tmp_path/'copy.cif';export_calculation_copy(structure,path)
    parser=CifParser(path,frac_tolerance=0);block=next(iter(parser.as_dict().values()))
    labels=block['_atom_site_label'];assert len(labels)==len(set(labels))
    reread=parser.parse_structures(primitive=False)[0]
    assert reread.composition.almost_equals(structure.composition)
    np.testing.assert_allclose(reread.frac_coords,structure.frac_coords,atol=1e-8)
    assert structure[0].label=='mixed'
