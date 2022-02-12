import numpy as np

from openmdao.api import n2
import openmdao.api as om

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.structures.struct_groups import SpatialBeamAlone
from openaerostruct.transfer.displacement_transfer_group import DisplacementTransferGroup

def StructualAnalysis(mesh_undeformed_left, wingInfo, loads):

    # Create the problem and assign the model group
    prob = om.Problem()

    # Add structural group to problem
    struct_group = SpatialBeamAlone(surface=wingInfo)
    prob.model.add_subsystem("struct_group", struct_group)

    # Add indep_vars to the structural group
    indep_var_comp = om.IndepVarComp()
    indep_var_comp.add_output("loads", val=loads, units="N")
    indep_var_comp.add_output("load_factor", val=1.0)
    struct_group.add_subsystem("indep_vars", indep_var_comp, promotes=["*"])

    # Add displacement transfer group
    disptrans_group = DisplacementTransferGroup(surface=wingInfo)
    prob.model.add_subsystem("disptransfer_group", disptrans_group)
    ny = wingInfo["mesh"].shape[1]
    disptrans_group.set_input_defaults('disp', val=np.ones((ny, 6)), units="m")

    prob.model.connect("struct_group.disp","disptransfer_group.disp")
    prob.model.connect("struct_group.mesh","disptransfer_group.mesh")
    prob.model.connect("struct_group.nodes","disptransfer_group.nodes")


    # Set up the problem
    prob.setup(force_alloc_complex=False)
    prob.run_model()
    # n2(prob)

    mesh_deformed_left = np.copy(prob["disptransfer_group.def_mesh"])

    # Decompose the mesh_deformed_left into x, y, z coordinates
    mesh_deformed_left_x = mesh_deformed_left[:,:,0]
    mesh_deformed_left_y = mesh_deformed_left[:,:,1]
    mesh_deformed_left_z = mesh_deformed_left[:,:,2]


    # mesh_delta_left = mesh_deformed_left - wingInfo["mesh"] # wingInfo["mesh"] is the initial mesh
    mesh_delta_left = mesh_deformed_left - mesh_undeformed_left
    mesh_delta_left_x = mesh_delta_left[:,:,0]
    mesh_delta_left_y = mesh_delta_left[:,:,1]
    mesh_delta_left_z = mesh_delta_left[:,:,2]


    # Compute the structural nodes based on the fem_origin location (weighted sum of the LE and TE mesh vertices)
    # Note that here the spar nodes are computed based on the mesh_deformed
    # In aerodynamics no structural analysis are conducted, these points are only used as locations to computeing the loads 
    # s_pts [ny, 3]
    # Taken from load_transfer.py
    # Note that the spar nodes depend ONLY on the LE and TE points
    # In fact all displacements along a rib is depend on the spar nodes 
    fem_origin = wingInfo["fem_origin"]
    spar_pts_deformed_left = (1 - fem_origin) * mesh_deformed_left[0, :, :] + fem_origin * mesh_deformed_left[-1, :, :]
    spar_pts_deformed_left_x = spar_pts_deformed_left[:,0]
    spar_pts_deformed_left_y = spar_pts_deformed_left[:,1]
    spar_pts_deformed_left_z = spar_pts_deformed_left[:,2]


    # Thickness
    # Tip to chord [m]
    sparThickness = np.copy(prob["struct_group.thickness"])

    # Radius
    # if Radius_cp not defined，radius wiil be dependent on t_over_c
    sparRadius = np.copy(prob["struct_group.radius"])

    # Intersection
    # If all the values should be negative, otherwie the hollow spar intersects itself 
    thickness_intersects = np.copy(prob["struct_group.thickness_intersects"])

    # Mass
    # Element mass
    element_mass = np.copy(prob["struct_group.element_mass"])

    # Total wing mass
    # Note that in the wingInfo there is a varaible call wing_weight_ratio
    wingStructuralMass = prob["struct_group.structural_mass"][0]	

    # Stess
    vonmisesStress = np.copy(prob["struct_group.vonmises"])

    # Failure
    # Negative -> Safe, Positive -> Fail
    failure = prob["struct_group.failure"][0]




    return(mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
        mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z, 
        spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
        sparThickness, sparRadius, thickness_intersects, 
        element_mass, wingStructuralMass, vonmisesStress, failure)