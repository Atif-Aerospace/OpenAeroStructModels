"""Test with using 
mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z
as inputs"""
import Aerodynamics.Aero_IntegratedWorkflow as aiw
import Aerodynamics.Aero_InputFileWimpress as aif # Read Input File

import Auxiliary.Aux_PlotWingAndLoads as xpw
import Auxiliary.Aux_AeroWriteToFile as xwf
import Auxiliary.Aux_ProduceThreeJS as xpt

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

mesh_delta_left_x = np.array([  [0.264464, 0.228736, 0.195599, 0.165053, 0.137098, 0.111735, 0.088963, 0.068783, 0.051194, 0.036198, 0.023163, 0.013025, 0.005785, 0.001443, -0.000000],
                                [0.264458, 0.228730, 0.195594, 0.165049, 0.137095, 0.111733, 0.088963, 0.068784, 0.051197, 0.036202, 0.023168, 0.013031, 0.005790, 0.001447, -0.000000],
                                [0.264453, 0.228725, 0.195589, 0.165044, 0.137092, 0.111732, 0.088963, 0.068786, 0.051200, 0.036206, 0.023173, 0.013036, 0.005795, 0.001450, -0.000000],
                                [0.264447, 0.228719, 0.195583, 0.165040, 0.137089, 0.111730, 0.088963, 0.068787, 0.051203, 0.036210, 0.023178, 0.013042, 0.005800, 0.001453, -0.000000],
                                [0.264441, 0.228713, 0.195578, 0.165036, 0.137086, 0.111728, 0.088963, 0.068788, 0.051206, 0.036214, 0.023183, 0.013047, 0.005805, 0.001456, -0.000000],
                                [0.264436, 0.228707, 0.195573, 0.165031, 0.137083, 0.111727, 0.088963, 0.068790, 0.051208, 0.036218, 0.023188, 0.013052, 0.005810, 0.001460, -0.000000],
                                [0.264430, 0.228702, 0.195567, 0.165027, 0.137080, 0.111725, 0.088962, 0.068791, 0.051211, 0.036222, 0.023193, 0.013058, 0.005815, 0.001463, -0.000000]])

mesh_delta_left_y = np.array([  [0.002645, 0.002287, 0.001956, 0.001651, 0.001371, 0.001117, 0.000890, 0.000688, 0.000512, 0.000362, 0.000232, 0.000130, 0.000058, 0.000014, -0.000000],
                                [0.002645, 0.002287, 0.001956, 0.001650, 0.001371, 0.001117, 0.000890, 0.000688, 0.000512, 0.000362, 0.000232, 0.000130, 0.000058, 0.000014, -0.000000],
                                [0.002645, 0.002287, 0.001956, 0.001650, 0.001371, 0.001117, 0.000890, 0.000688, 0.000512, 0.000362, 0.000232, 0.000130, 0.000058, 0.000014, -0.000000],
                                [0.002644, 0.002287, 0.001956, 0.001650, 0.001371, 0.001117, 0.000890, 0.000688, 0.000512, 0.000362, 0.000232, 0.000130, 0.000058, 0.000015, -0.000000],
                                [0.002644, 0.002287, 0.001956, 0.001650, 0.001371, 0.001117, 0.000890, 0.000688, 0.000512, 0.000362, 0.000232, 0.000130, 0.000058, 0.000015, -0.000000],
                                [0.002644, 0.002287, 0.001956, 0.001650, 0.001371, 0.001117, 0.000890, 0.000688, 0.000512, 0.000362, 0.000232, 0.000131, 0.000058, 0.000015, -0.000000],
                                [0.002644, 0.002287, 0.001956, 0.001650, 0.001371, 0.001117, 0.000890, 0.000688, 0.000512, 0.000362, 0.000232, 0.000131, 0.000058, 0.000015, -0.000000]])

mesh_delta_left_z = np.array([  [2.644637, 2.287360, 1.955991, 1.650531, 1.370983, 1.117348, 0.889629, 0.687826, 0.511943, 0.361980, 0.231631, 0.130253, 0.057852, 0.014433, -0.000000],
                                [2.644581, 2.287303, 1.955938, 1.650487, 1.370952, 1.117332, 0.889628, 0.687841, 0.511971, 0.362020, 0.231680, 0.130307, 0.057902, 0.014465, -0.000000],
                                [2.644525, 2.287245, 1.955885, 1.650444, 1.370921, 1.117316, 0.889627, 0.687856, 0.512000, 0.362059, 0.231730, 0.130361, 0.057951, 0.014498, -0.000000],
                                [2.644469, 2.287188, 1.955832, 1.650400, 1.370890, 1.117299, 0.889627, 0.687870, 0.512028, 0.362098, 0.231779, 0.130415, 0.058000, 0.014531, -0.000000],
                                [2.644413, 2.287131, 1.955779, 1.650356, 1.370859, 1.117283, 0.889626, 0.687885, 0.512056, 0.362137, 0.231829, 0.130469, 0.058050, 0.014563, -0.000000],
                                [2.644357, 2.287073, 1.955727, 1.650313, 1.370828, 1.117266, 0.889625, 0.687899, 0.512085, 0.362177, 0.231878, 0.130523, 0.058099, 0.014596, -0.000000],
                                [2.644301, 2.287016, 1.955674, 1.650269, 1.370796, 1.117250, 0.889624, 0.687914, 0.512113, 0.362216, 0.231928, 0.130577, 0.058149, 0.014628, -0.000000]])

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
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
    ) = aiw.AerodynamicsFunction_wp_10_DM(aif.wing_area_wpref, aif.aspect_ratio_wpref, aif.kink_location_ratio, aif.body_side_ratio,   # Planform definition  
    aif.taper_ratio_trap, aif.root_chord_extension_ratio, aif.trap_quarter_sweep,   # Planform definition
    aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
    aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
    aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
    mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z)    # these has to be assumed

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

produceThreeJS = True
if produceThreeJS == True:
    fileName = "ThreeJsPlot\models\OnlyAero\OnlyAeroDesignPoint_AoA_" + f"{aif.alpha:.1f}" + "_Mach_" + f"{aif.Mach_number:.1f}"+ ".json"
    xpt.ThreeJSFile(fileName, mesh_deformed_left, mesh_point_force, mesh_delta_left)


# write to file
fileName = "Results/Model_Execution/AerodynamicsDesignReport_Wim10_DM.txt"
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


# Produce Drag Polar
dragPolar = True

if dragPolar == True:

    print("Producing Drag Polar\n")

    # Create arrays to hold polar data
    N_Mach = 10
    N_alpha = 16
    N_CL = 10

    Mach_dragPolar = np.linspace(0.1, 1.0, N_Mach)
    alpha_dragPolar = np.linspace(-5, 10, N_alpha)

    CL_dragPolar_Raw = np.zeros((N_Mach, N_alpha))
    CD_total_dragPolar_Raw = np.zeros((N_Mach, N_alpha))
    CD_i_dragPolar_raw = np.zeros((N_Mach, N_alpha))

    # CL_dragPolar_fixed = np.linspace(0.0, 0.50, N_CL)
    CL_dragPolar_fixed = np.array([0.00000, 0.05000, 0.10000, 0.15000, 0.20000, 0.25000, 0.30000, 0.40000, 0.50000, 0.60000]) 
    CD_dragPolar_intp = np.zeros((N_Mach, N_CL))
    CD_i_dragPolar_intp = np.zeros((N_Mach, N_CL))

    # Sweep through alphas and Mach to create polar
    for i in range(len(Mach_dragPolar)):

        for j in range(len(alpha_dragPolar)):

            # Set Mach
            aif.Mach_number = Mach_dragPolar[i] 

            # Set alpha
            aif.alpha = alpha_dragPolar[j]

            # These need more checks   
            aif.velocity = Mach_dragPolar[i] * 296.54 # scale the speed
            aif.Re = aif.velocity / 248.136000 * 1000000 # scale the Reynolds number
            
            print("Mach = ", aif.Mach_number, ", alpha = ", aif.alpha, "\n")
            # Run analysis
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
            mesh_delta_left_dragPolar, # mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
            mesh_deformed_left_dragPolar, mesh_deformed_left_x_dragPolar, mesh_deformed_left_y_dragPolar, mesh_deformed_left_z_dragPolar,
            spar_pts_deformed_left_dragPolar, spar_pts_deformed_left_x_dragPolar, spar_pts_deformed_left_y_dragPolar, spar_pts_deformed_left_z_dragPolar,
            wing_area_asref_dragPolar, L_Wing_dragPolar, D_Wing_dragPolar, LoD_Wing_dragPolar,
            CL_Wing_total_dragPolar, CD_Wing_i_dragPolar, CD_Wing_v_dragPolar, CD_Wing_w_dragPolar, CD_Wing_total_dragPolar, 
            CM_Wing_dragPolar, CM_Wing_roll_dragPolar, CM_Wing_pitch_dragPolar, CM_Wing_yaw_dragPolar,
            sec_forces_dragPolar, sec_forces_Fx_dragPolar, sec_forces_Fy_dragPolar, sec_forces_Fz_dragPolar, mesh_point_force_dragPolar,
            loads_dragPolar, loads_Fx_dragPolar, loads_Fy_dragPolar, loads_Fz_dragPolar, loads_Mx_dragPolar, loads_My_dragPolar, loads_Mz_dragPolar
            ) = aiw.AerodynamicsFunction_wp_10_DM(aif.wing_area_wpref, aif.aspect_ratio_wpref, aif.kink_location_ratio, aif.body_side_ratio,   # Planform definition  
            aif.taper_ratio_trap, aif.root_chord_extension_ratio, aif.trap_quarter_sweep,   # Planform definition
            aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
            aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
            aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
            aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
            mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z)

            # Record CL, CD
            CL_dragPolar_Raw[i,j] = CL_Wing_total_dragPolar
            CD_total_dragPolar_Raw[i,j] = CD_Wing_total_dragPolar
            CD_i_dragPolar_raw[i,j] = CD_Wing_i_dragPolar # induced drag
                        

            produceThreeJS = True
            if produceThreeJS == True:
                fileName = "ThreeJsPlot\models\OnlyAero\DragPolar\OnlyAeroDragPolar_AoA_" +  f"{alpha_dragPolar[j]:.1f}" + "_Mach_" + f"{Mach_dragPolar[i]:.1f}" + ".json"
                xpt.ThreeJSFile(fileName, mesh_deformed_left_dragPolar, mesh_point_force_dragPolar, mesh_delta_left_dragPolar)
                # Note that the mesh_deformed_left, and mesh_delta_left should actually change with the aoa and mach
                # they are currently fixed

        # Using Interpolate to get CD at fixed CLs
        tck_total = interpolate.splrep(CL_dragPolar_Raw[i,:], CD_total_dragPolar_Raw[i,:], s=0)
        CD_dragPolar_intp[i,:] = interpolate.splev(CL_dragPolar_fixed, tck_total, der=0)

        tck_induced = interpolate.splrep(CL_dragPolar_Raw[i,:], CD_i_dragPolar_raw[i,:], s=0)
        CD_i_dragPolar_intp[i,:] = interpolate.splev(CL_dragPolar_fixed, tck_induced, der=0)

    # Plot polar
    plotSwitch = True
    if plotSwitch == True:
        plt.figure()

        # Sweep through Mach to create polar
        for i in range(len(Mach_dragPolar)):

            plt.plot(CD_total_dragPolar_Raw[i,:] * 1e4, CL_dragPolar_Raw[i,:], 'x')
            plt.plot(CD_dragPolar_intp[i,:] * 1e4, CL_dragPolar_fixed, "-o")
            plt.grid(color="lightgray", linestyle="-", linewidth=1)

        # plt.legend(['True', 'Interpolate', 'True'])
        plt.xlabel("$C_D$ (counts)")
        plt.ylabel("$C_L$")
        plt.title("Drag polar")
        plt.show()
        # plt.show(block=False)