import numpy as np

import Aerodynamics.Aero_IntegratedWorkflow as aiw
import Aerodynamics.Aero_InputFileWimpress as aif

import Structures.Struct_IntegratedWorkflow as siw
import Structures.Struct_InputFileWimpress as sif

import Auxiliary.Aux_PlotWingAndLoads as xpw
import Auxiliary.Aux_AeroWriteToFile as xwf

def AddNumbers(x1, x2):
	return (x1 + x2)

def MinusNumbers(x1, x2):
	return (x1 - x2)

def MultiplyNumbers(x1, x2):
	return (x1 * x2)

def TwoOutputsModel(x1, x2):
	y1 = x1 + x2
	y2 = x1 * x2
	return (y1, y2)
	
def OAS_Aerodynamics_V1(SW, AR, Kink, TR, Sweep, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z):

    aif.wing_area_wpref = SW
    aif.aspect_ratio_wpref = AR
    aif.kink_location_ratio = Kink
    aif.taper_ratio_trap = TR
    aif.trap_quarter_sweep = Sweep
    
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
        mesh_delta_left, # mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
        mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
        spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
        wing_area_asref, L_Wing, D_Wing, LoD_Wing,
        CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
        CL_DP_fixed, CD_DP_final, CD_i_DP_final,
        CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
        sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
        loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
        ) = aiw.AerodynamicsFunction_wp_10_DM(aif.wing_area_wpref, aif.aspect_ratio_wpref, aif.kink_location_ratio, aif.body_side_ratio,   # Planform definition  
        aif.taper_ratio_trap, aif.root_chord_extension_ratio, aif.trap_quarter_sweep,   # Planform definition
        aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
        aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
        aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
        aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
        mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z)

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
    plotSwitch = False
    if plotSwitch == True:
        xpw.PlotWingAndLoads(mesh_initial_left, mesh_undeformed_left, mesh_deformed_left,
        spar_pts_undeformed_left, spar_pts_deformed_left,
        loads)

    # write to file
    # fileName = "Results/Model_Execution/AerodynamicsDesignReport_Wim10_DM.txt"
    # xwf.WriteToFiles(fileName, aif.wing_area_wpref, aif.aspect_ratio_wpref, wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap,  # Planform definition 
    #     aif.kink_location_ratio, aif.body_side_ratio, taper_ratio_t_r, taper_ratio_k_r, aif.taper_ratio_trap, aif.root_chord_extension_ratio,    # Planform definition
    #     entire_span, half_span, kink_location, root_chord_trap, root_chord, kink_chord, tip_chord, 
    #     inboard_quarter_sweep, outboard_quarter_sweep, inboard_LE_sweep, outboard_LE_sweep, aif.trap_quarter_sweep,
    #     aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    #     aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
    #     aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
    #     aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
    #     aif.delta_x_LE_1, aif.delta_x_LE_2, aif.delta_x_LE_3,       # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
    #     aif.delta_y_LE_1, aif.delta_y_LE_2, aif.delta_y_LE_3,
    #     aif.delta_z_LE_1, aif.delta_z_LE_2, aif.delta_z_LE_3,
    #     aif.delta_x_TE_1, aif.delta_x_TE_2, aif.delta_x_TE_3,       # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
    #     aif.delta_y_TE_1, aif.delta_y_TE_2, aif.delta_y_TE_3,
    #     aif.delta_z_TE_1, aif.delta_z_TE_2, aif.delta_z_TE_3,
    #     wingarea_wpref_dependency, aspectratio_wpref_dependency, 
    #     wingarea_geo_dependency, aspectratio_geo_dependency, 
    #     wingarea_trap_dependency, aspectratio_trap_dependency,
    #     kinklocation_ratio_dependency, bodysideratio_dependency, rootchord_extratio_dependency, taperratioTrap_dependency, taperratioTr_dependency, taperratioKr_dependency,
    #     entirespan_dependency, halfspan_dependency, kinklocation_dependency, rootchordTrap_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency,
    #     trapquartersweep_dependency, inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency,
    #     nx, ny_inboard, ny_outboard, ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z,
    #     t_over_c_actual, twist_actual, wingInfo,
    #     mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    #     spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z,
    #     mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    #     mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    #     spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
    #     wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    #     CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
    #     CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    #     sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
    #     loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz)

    return (LoD_Wing, CL_Wing_total, CL_DP_fixed, CD_DP_final, CD_i_DP_final, 
    loads, # If we need to plot uncomment the lines below
    mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z,
    spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z
    )

def OAS_Structure_V1(SW, AR, Kink, TR, Sweep, loads):

    aif.wing_area_wpref = SW
    aif.aspect_ratio_wpref = AR
    aif.kink_location_ratio = Kink
    aif.taper_ratio_trap = TR
    aif.trap_quarter_sweep = Sweep
    sif.loads = loads

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
    t_over_c_actual, twist_actual, wingInfo, loads,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z,
    mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
    sparThickness, sparRadius, thickness_intersects, 
    element_mass, wingStructuralMass, vonmisesStress, failure
    ) = siw.StructFunction_wp_10(sif.E, sif.G, sif.yieldStress, sif.mrho, sif.fem_origin, sif.wing_weight_ratio, # Assumed variables
    sif.wing_area_wpref, sif.aspect_ratio_wpref, sif.kink_location_ratio, sif.body_side_ratio, sif.taper_ratio_trap, sif.root_chord_extension_ratio, sif.trap_quarter_sweep,   # Planform definition
    sif.dihedral, sif.twist_cp_1, sif.twist_cp_2, sif.twist_cp_3, sif.t_over_c_cp_1, sif.t_over_c_cp_2, sif.t_over_c_cp_3, sif.c_max_t,  # 3D geometry: 1, 2, 3 are from tip to root
    sif.thickness_cp_1, sif.thickness_cp_2, sif.thickness_cp_3, #thickness control point
    sif.loads # sif.loads_Fx, sif.loads_Fy, sif.loads_Fz, sif.loads_Mx, sif.loads_My, sif.loads_Mz # Loads
    )

    # Plot
    plotSwitch = False
    if plotSwitch == True:
        xpw.PlotWingAndLoads(mesh_initial_left, mesh_undeformed_left, mesh_deformed_left,
        spar_pts_undeformed_left, spar_pts_deformed_left,
        loads)
    
    return (thickness_intersects, element_mass, wingStructuralMass, vonmisesStress, failure,
    # mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z)     # If we need to plot comment this line and uncomment the lines below
    mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z, 
    spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z,
    spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z
    )

def OAS_Iteration_V1(SW, AR, Kink, TR, Sweep):

    aif.wing_area_wpref = SW
    aif.aspect_ratio_wpref = AR
    aif.kink_location_ratio = Kink
    aif.taper_ratio_trap = TR
    aif.trap_quarter_sweep = Sweep

    # Intial assumed deformation for Aero
    mesh_delta_left_x_Aero = np.zeros((7, 15))
    mesh_delta_left_y_Aero = np.zeros((7, 15))
    mesh_delta_left_z_Aero = np.zeros((7, 15))


    # Intial assumed deformation for Struct
    mesh_delta_left_x_Struct = np.zeros((7, 15))
    mesh_delta_left_y_Struct = np.zeros((7, 15))
    mesh_delta_left_z_Struct = np.zeros((7, 15))

    # Count number of iterations
    iterCount = 0

    # maximum iteration
    iterMax = 100

    # gap between the deformed mesh from aero and structure, initial error is set to 1 to start the while loop
    error = 1

    # convergence tolerance
    tol = 1e-6

    # iteration to reduce the gap
    while (error > tol and iterCount < iterMax):

        # pass deformation from Struct to Aero
        mesh_delta_left_x_Aero = np.copy(mesh_delta_left_x_Struct)
        mesh_delta_left_y_Aero = np.copy(mesh_delta_left_y_Struct)
        mesh_delta_left_z_Aero = np.copy(mesh_delta_left_z_Struct)

        # Run Aero
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
            mesh_undeformed_left_Aero, mesh_undeformed_left_x_Aero, mesh_undeformed_left_y_Aero, mesh_undeformed_left_z_Aero,
            spar_pts_undeformed_left_Aero, spar_pts_undeformed_left_x_Aero, spar_pts_undeformed_left_y_Aero, spar_pts_undeformed_left_z_Aero,
            mesh_delta_left_Aero, # mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
            mesh_deformed_left_Aero, mesh_deformed_left_x_Aero, mesh_deformed_left_y_Aero, mesh_deformed_left_z_Aero,
            spar_pts_deformed_left_Aero, spar_pts_deformed_left_x_Aero, spar_pts_deformed_left_y_Aero, spar_pts_deformed_left_z_Aero,
            wing_area_asref, L_Wing, D_Wing, LoD_Wing,
            CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
            CL_DP_fixed, CD_DP_final, CD_i_DP_final,
            CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
            sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
            loads_Aero, loads_Fx_Aero, loads_Fy_Aero, loads_Fz_Aero, loads_Mx_Aero, loads_My_Aero, loads_Mz_Aero
            ) = aiw.AerodynamicsFunction_wp_10_DM(aif.wing_area_wpref, aif.aspect_ratio_wpref, aif.kink_location_ratio, aif.body_side_ratio,   # Planform definition  
            aif.taper_ratio_trap, aif.root_chord_extension_ratio, aif.trap_quarter_sweep,   # Planform definition
            aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
            aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
            aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
            aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
            mesh_delta_left_x_Aero, mesh_delta_left_y_Aero, mesh_delta_left_z_Aero)

        # Pass loads from Aero to Struct
        loads_Struct = np.copy(loads_Aero)

        # Run Struct
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
            t_over_c_actual, twist_actual, wingInfo, loads_Struct,
            mesh_undeformed_left_Struct, mesh_undeformed_left_x_Struct, mesh_undeformed_left_y_Struct, mesh_undeformed_left_z_Struct,
            spar_pts_undeformed_left_Struct, spar_pts_undeformed_left_x_Struct, spar_pts_undeformed_left_y_Struct, spar_pts_undeformed_left_z_Struct,
            mesh_delta_left_Struct, mesh_delta_left_x_Struct, mesh_delta_left_y_Struct, mesh_delta_left_z_Struct,
            mesh_deformed_left_Struct, mesh_deformed_left_x_Struct, mesh_deformed_left_y_Struct, mesh_deformed_left_z_Struct,
            spar_pts_deformed_left_Struct, spar_pts_deformed_left_x_Struct, spar_pts_deformed_left_y_Struct, spar_pts_deformed_left_z_Struct,
            sparThickness, sparRadius, thickness_intersects, 
            element_mass, wingStructuralMass, vonmisesStress, failure
            ) = siw.StructFunction_wp_10(sif.E, sif.G, sif.yieldStress, sif.mrho, sif.fem_origin, sif.wing_weight_ratio, # Assumed variables
            sif.wing_area_wpref, sif.aspect_ratio_wpref, sif.kink_location_ratio, sif.body_side_ratio, sif.taper_ratio_trap, sif.root_chord_extension_ratio, sif.trap_quarter_sweep,   # Planform definition
            sif.dihedral, sif.twist_cp_1, sif.twist_cp_2, sif.twist_cp_3, sif.t_over_c_cp_1, sif.t_over_c_cp_2, sif.t_over_c_cp_3, sif.c_max_t,  # 3D geometry: 1, 2, 3 are from tip to root
            sif.thickness_cp_1, sif.thickness_cp_2, sif.thickness_cp_3, #thickness control point
            loads_Struct # sif.loads_Fx, sif.loads_Fy, sif.loads_Fz, sif.loads_Mx, sif.loads_My, sif.loads_Mz # Loads
            )

        # current gap between the deformed mesh from aero and structure
        error = np.square(np.subtract(mesh_deformed_left_Struct, mesh_deformed_left_Aero)).mean()

        # Count the number of iterations
        iterCount += 1

        print(iterCount,": ",error)
        

    # Final converged values
    loads_Final = np.copy(loads_Aero)
    loads_Fx_Final = np.copy(loads_Fx_Aero)
    loads_Fy_Final = np.copy(loads_Fy_Aero)
    loads_Fz_Final = np.copy(loads_Fz_Aero)
    loads_Mx_Final = np.copy(loads_Mx_Aero)
    loads_My_Final = np.copy(loads_My_Aero)
    loads_Mz_Final = np.copy(loads_Mz_Aero)

    mesh_undeformed_left_Final = np.copy(mesh_undeformed_left_Struct)
    mesh_undeformed_left_x_Final = np.copy(mesh_undeformed_left_x_Struct)
    mesh_undeformed_left_y_Final = np.copy(mesh_undeformed_left_y_Struct)
    mesh_undeformed_left_z_Final = np.copy(mesh_undeformed_left_z_Struct)

    mesh_deformed_left_Final = np.copy(mesh_deformed_left_Struct)
    mesh_deformed_left_x_Final = np.copy(mesh_deformed_left_x_Struct)
    mesh_deformed_left_y_Final = np.copy(mesh_deformed_left_y_Struct)
    mesh_deformed_left_z_Final = np.copy(mesh_deformed_left_z_Struct)

    mesh_delta_left_Final = np.copy(mesh_delta_left_Struct)
    mesh_delta_left_x_Final = np.copy(mesh_delta_left_x_Struct)
    mesh_delta_left_y_Final = np.copy(mesh_delta_left_y_Struct)
    mesh_delta_left_z_Final = np.copy(mesh_delta_left_z_Struct)

    spar_pts_undeformed_left_Final = np.copy(spar_pts_undeformed_left_Struct)
    spar_pts_undeformed_left_x_Final = np.copy(spar_pts_undeformed_left_x_Struct)
    spar_pts_undeformed_left_y_Final = np.copy(spar_pts_undeformed_left_y_Struct)
    spar_pts_undeformed_left_z_Final = np.copy(spar_pts_undeformed_left_z_Struct)

    spar_pts_deformed_left_Final = np.copy(spar_pts_deformed_left_Struct)
    spar_pts_deformed_left_x_Final = np.copy(spar_pts_deformed_left_x_Struct)
    spar_pts_deformed_left_y_Final = np.copy(spar_pts_deformed_left_y_Struct) 
    spar_pts_deformed_left_z_Final = np.copy(spar_pts_deformed_left_z_Struct)

    # Print
    printSwitch = False
    if printSwitch == True:
        print("================================= Results ===================================\n") 
        print("Final aerodynamics loads:\n", loads_Final,"\n") 
        print("Final deformed mesh:\n", mesh_deformed_left_Final,"\n")  
        print("aeropoint_group.aero_states.wing_sec_forces:\n", sec_forces,"\n") 
        print("loadtransfer_group.sec_forces:\n", sec_forces,"\n") 

    # Temp
    # import numpy as np
    # np.savetxt("mesh_deformed_left_x.csv", mesh_deformed_left_x, delimiter=",")
    # np.savetxt("mesh_deformed_left_y.csv", mesh_deformed_left_y, delimiter=",")
    # np.savetxt("mesh_deformed_left_z.csv", mesh_deformed_left_z, delimiter=",")

    # Plot
    plotSwitch = False
    if plotSwitch == True:
        xpw.PlotWingAndLoads(mesh_initial_left, mesh_undeformed_left_Final, mesh_deformed_left_Final,
        spar_pts_undeformed_left_Final, spar_pts_deformed_left_Final,
        loads_Final)

    # write to file
    fileName = "Results/Model_Execution/CombinedDesignReport_Iteration_FineMesh.txt"
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
        mesh_undeformed_left_Final, mesh_undeformed_left_x_Final, mesh_undeformed_left_y_Final, mesh_undeformed_left_z_Final,
        spar_pts_undeformed_left_Final, spar_pts_undeformed_left_x_Final, spar_pts_undeformed_left_y_Final, spar_pts_undeformed_left_z_Final,
        mesh_delta_left_Final, mesh_delta_left_x_Final, mesh_delta_left_y_Final, mesh_delta_left_z_Final,
        mesh_deformed_left_Final, mesh_deformed_left_x_Final, mesh_deformed_left_y_Final, mesh_deformed_left_z_Final,
        spar_pts_deformed_left_Final, spar_pts_deformed_left_x_Final, spar_pts_deformed_left_y_Final, spar_pts_deformed_left_z_Final,
        wing_area_asref, L_Wing, D_Wing, LoD_Wing,
        CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
        CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
        sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
        loads_Final, loads_Fx_Final, loads_Fy_Final, loads_Fz_Final, loads_Mx_Final, loads_My_Final, loads_Mz_Final)

    return (LoD_Wing, CL_Wing_total, CL_DP_fixed, CD_DP_final, CD_i_DP_final, 
    thickness_intersects, element_mass, wingStructuralMass, vonmisesStress, failure, # If we need to plot uncomment the lines below
    loads_Final, 
    mesh_undeformed_left_x_Final, mesh_undeformed_left_y_Final, mesh_undeformed_left_z_Final,
    mesh_deformed_left_x_Final, mesh_deformed_left_y_Final, mesh_deformed_left_z_Final,
    mesh_delta_left_x_Final, mesh_delta_left_y_Final, mesh_delta_left_z_Final, 
    spar_pts_undeformed_left_x_Final, spar_pts_undeformed_left_y_Final, spar_pts_undeformed_left_z_Final,
    spar_pts_deformed_left_x_Final, spar_pts_deformed_left_y_Final, spar_pts_deformed_left_z_Final)