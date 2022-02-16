"""this input file is used for validation, the values are changed to match the example from 
http://mdolab.engin.umich.edu/OpenAeroStruct/aero_walkthrough.html"""

import Aerodynamics.Aero_IntegratedWorkflow as aiw
import Aerodynamics.Aero_InputFileWimpress_Validation as aif # Read Input File
import Auxiliary.Aux_PlotWingAndLoads as xpw
import Auxiliary.Aux_AeroWriteToFile as xwf

# Call the Integrated Workflow as a function
(wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap, entire_span, half_span, kink_location,
    root_chord_trap, root_chord, kink_chord, tip_chord, taper_ratio_t_r, taper_ratio_k_r,
    inboard_quarter_sweep, outboard_quarter_sweep, inboard_LE_sweep, outboard_LE_sweep,
    wingarea_wpref_dependency, aspectratio_wpref_dependency, 
    wingarea_geo_dependency, aspectratio_geo_dependency, 
    wingarea_trap_dependency, aspectratio_trap_dependency,
    kinklocation_ratio_dependency, bodysideratio_dependency, rootchord_extratio_dependency, taperratioTrap_dependency, taperratioTr_dependency, taperratioKr_dependency,
    entirespan_dependency, halfspan_dependency, kinklocation_dependency, rootchordTrap_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency,
    trapquartersweep_dependency, inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency,
    nx, ny_inboard, ny_outboard, ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z,
    t_over_c_actual, twist_actual, wingInfo,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z,
    mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
    wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
    ) = aiw.AerodynamicsFunction_wp_10(aif.wing_area_wpref, aif.aspect_ratio_wpref, aif.kink_location_ratio, aif.body_side_ratio,   # Planform definition  
    aif.taper_ratio_trap, aif.root_chord_extension_ratio, aif.trap_quarter_sweep,   # Planform definition
    aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
    aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
    aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
    aif.delta_x_LE_1, aif.delta_x_LE_2, aif.delta_x_LE_3,   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
    aif.delta_y_LE_1, aif.delta_y_LE_2, aif.delta_y_LE_3,
    aif.delta_z_LE_1, aif.delta_z_LE_2, aif.delta_z_LE_3,
    aif.delta_x_TE_1, aif.delta_x_TE_2, aif.delta_x_TE_3,   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
    aif.delta_y_TE_1, aif.delta_y_TE_2, aif.delta_y_TE_3,
    aif.delta_z_TE_1, aif.delta_z_TE_2, aif.delta_z_TE_3)

# Print
printSwitch = False
if printSwitch == True:
    print("================================= Results ===================================\n") 
    print("Final aerodynamics loads:\n", loads,"\n") 
    print("Final deformed mesh:\n", mesh_deformed_left,"\n")  
    print("aeropoint_group.aero_states.wing_sec_forces:\n", sec_forces,"\n") 
    print("loadtransfer_group.sec_forces:\n", sec_forces,"\n") 

# Temp
# import numpy as np
# np.savetxt("mesh_deformed_left_x.csv", mesh_deformed_left_x, delimiter=",")
# np.savetxt("mesh_deformed_left_y.csv", mesh_deformed_left_y, delimiter=",")
# np.savetxt("mesh_deformed_left_z.csv", mesh_deformed_left_z, delimiter=",")

# Plot
plotSwitch = True
if plotSwitch == True:
    xpw.PlotWingAndLoads(mesh_initial_left, mesh_undeformed_left, mesh_deformed_left,
    spar_pts_undeformed_left, spar_pts_deformed_left,
    loads)

# write to file
fileName = "Results/Model_Execution/AerodynamicsDesignReport_Wim10_Validation.txt"
xwf.WriteToFiles(fileName, aif.wing_area_wpref, aif.aspect_ratio_wpref, wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap,  # Planform definition 
    aif.kink_location_ratio, aif.body_side_ratio, taper_ratio_t_r, taper_ratio_k_r, aif.taper_ratio_trap, aif.root_chord_extension_ratio,    # Planform definition
    entire_span, half_span, kink_location, root_chord_trap, root_chord, kink_chord, tip_chord, 
    inboard_quarter_sweep, outboard_quarter_sweep, inboard_LE_sweep, outboard_LE_sweep, aif.trap_quarter_sweep,
    aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
    aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
    aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
    aif.delta_x_LE_1, aif.delta_x_LE_2, aif.delta_x_LE_3,       # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
    aif.delta_y_LE_1, aif.delta_y_LE_2, aif.delta_y_LE_3,
    aif.delta_z_LE_1, aif.delta_z_LE_2, aif.delta_z_LE_3,
    aif.delta_x_TE_1, aif.delta_x_TE_2, aif.delta_x_TE_3,       # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
    aif.delta_y_TE_1, aif.delta_y_TE_2, aif.delta_y_TE_3,
    aif.delta_z_TE_1, aif.delta_z_TE_2, aif.delta_z_TE_3,
    wingarea_wpref_dependency, aspectratio_wpref_dependency, 
    wingarea_geo_dependency, aspectratio_geo_dependency, 
    wingarea_trap_dependency, aspectratio_trap_dependency,
    kinklocation_ratio_dependency, bodysideratio_dependency, rootchord_extratio_dependency, taperratioTrap_dependency, taperratioTr_dependency, taperratioKr_dependency,
    entirespan_dependency, halfspan_dependency, kinklocation_dependency, rootchordTrap_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency,
    trapquartersweep_dependency, inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency,
    nx, ny_inboard, ny_outboard, ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z,
    t_over_c_actual, twist_actual, wingInfo,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z,
    mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
    wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz)
