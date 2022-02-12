"""
Produce underformed mesh using OpenAeroStruct Geometry Component
This mesh has twist and dihedral, but no structural deformation
"""

from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint
from openaerostruct.structures.wingbox_fuel_vol_delta import WingboxFuelVolDelta
import openmdao.api as om

def UndeformedMesh(mesh_initial_right, dihedral, twist_cp, t_over_c_cp, CL0, CD0, k_lam, c_max_t, fem_origin):
    wingInfo = {
        # Wing definition
        "name": "wing",  # name of the surface
        "symmetry": True,  # if true, model one half of wing reflected across the plane y = 0
        "S_ref_type": "wetted",  # how we compute the wing area, can be 'wetted' or 'projected'
        "mesh": mesh_initial_right,
        "fem_model_type": "tube",
        "dihedral": dihedral,
        "twist_cp": twist_cp,
        # "thickness_cp": thickness_cp,
        "t_over_c_cp": t_over_c_cp, 
        "CL0": CL0,  # CL of the surface at alpha=0
        "CD0": CD0,  # CD of the surface at alpha=0
        # Airfoil properties for viscous drag calculation
        "k_lam": k_lam,  # percentage of chord with laminar flow, used for viscous drag 
        "c_max_t": c_max_t,  # chordwise location of maximum thickness
        "with_viscous": True,
        "with_wave": True,  # if true, compute wave drag
        # # Structural values are based on aluminum 7075
        # "E": E,  # [Pa] Young's modulus of the spar
        # "G": G,  # [Pa] shear modulus of the spar
        # "yield": yieldStress,  # [Pa] yield stress divided by 2.5 for limiting case
        # "mrho": mrho,  # [kg/m^3] material density
        "fem_origin": fem_origin,  # normalized chordwise location of the spar
        # "wing_weight_ratio": wing_weight_ratio,
        "struct_weight_relief": False,  # True to add the weight of the structure to the loads on the structure
        "distributed_fuel_weight": False,
        # Constraints
        "exact_failure_constraint": False,  # if false, use KS function
    }


    # Create the OpenMDAO problem
    # This prob is used to produce the mesh with twist (no structural deformation at the moment)
    prob_mesh = om.Problem()

    # Create and add a group that handles the geometry for the
    # aerodynamic lifting surface
    winggeom_group = Geometry(surface=wingInfo)
    prob_mesh.model.add_subsystem("winggeom_group", winggeom_group)

    # Set up and run the model
    prob_mesh.setup()
    prob_mesh.run_model()
    # n2(prob_mesh)

    # Thickness to chord (from tip to root):
    t_over_c_actual = prob_mesh["winggeom_group.t_over_c"]
    # Actual_twist [deg] (from tip to root):
    twist_actual = prob_mesh["winggeom_group.twist"]

    # mesh after twist and dihedral?, no structural deformation in aerodynamics analysis)
    mesh_undeformed_left = prob_mesh["winggeom_group.mesh"]
    mesh_undeformed_left_x = mesh_undeformed_left[:,:,0]
    mesh_undeformed_left_y = mesh_undeformed_left[:,:,1]
    mesh_undeformed_left_z = mesh_undeformed_left[:,:,2]

    spar_pts_undeformed_left = (1 - fem_origin) * mesh_undeformed_left[0, :, :] + fem_origin * mesh_undeformed_left[-1, :, :]
    spar_pts_undeformed_left_x = spar_pts_undeformed_left[:,0]
    spar_pts_undeformed_left_y = spar_pts_undeformed_left[:,1]
    spar_pts_undeformed_left_z = spar_pts_undeformed_left[:,2]

    return (t_over_c_actual, twist_actual, wingInfo,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z)