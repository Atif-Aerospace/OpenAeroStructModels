"""A combined single run with aero and struct loop solved by iterations
note that we need to change both aif and sif to make the palnform value the same
sif.wing_area_wpref, sif.aspect_ratio_wpref may not be equal to aif.wing_area_wpref, aif.aspect_ratio_wpref"""

import Aerodynamics.Aero_IntegratedWorkflow as aiw
import Aerodynamics.Aero_InputFileWimpress as aif # Read Input File

import Structures.Struct_IntegratedWorkflow as siw
import Structures.Struct_InputFileWimpress as sif

import Auxiliary.Aux_PlotWingAndLoads as xpw
import Auxiliary.Aux_AeroWriteToFile as xwf
import Auxiliary.Aux_ProduceThreeJS as xpt

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

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

print("Design Point")
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
        CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
        sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
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
        ) = siw.StructFunction_wp_10(
        sif.wing_area_wpref, sif.aspect_ratio_wpref, sif.kink_location_ratio, sif.body_side_ratio, sif.taper_ratio_trap, sif.root_chord_extension_ratio, sif.trap_quarter_sweep,   # Planform definition
        sif.dihedral, sif.twist_cp_1, sif.twist_cp_2, sif.twist_cp_3, sif.t_over_c_cp_1, sif.t_over_c_cp_2, sif.t_over_c_cp_3, sif.c_max_t,  # 3D geometry: 1, 2, 3 are from tip to root
        sif.thickness_cp_1, sif.thickness_cp_2, sif.thickness_cp_3, #thickness control point
        sif.E, sif.G, sif.yieldStress, sif.mrho, sif.fem_origin, sif.wing_weight_ratio, # Assumed variables
        loads_Struct # sif.loads_Fx, sif.loads_Fy, sif.loads_Fz, sif.loads_Mx, sif.loads_My, sif.loads_Mz # Loads
        )

    # current gap between the deformed mesh from aero and structure
    error = np.square(np.subtract(mesh_deformed_left_Struct, mesh_deformed_left_Aero)).mean()

    # Count the number of iterations
    iterCount += 1

    print(iterCount,": ",error)
    

# Converged converged values
loads_Converged = np.copy(loads_Aero)
loads_Fx_Converged = np.copy(loads_Fx_Aero)
loads_Fy_Converged = np.copy(loads_Fy_Aero)
loads_Fz_Converged = np.copy(loads_Fz_Aero)
loads_Mx_Converged = np.copy(loads_Mx_Aero)
loads_My_Converged = np.copy(loads_My_Aero)
loads_Mz_Converged = np.copy(loads_Mz_Aero)

mesh_undeformed_left_Converged = np.copy(mesh_undeformed_left_Struct)
mesh_undeformed_left_x_Converged = np.copy(mesh_undeformed_left_x_Struct)
mesh_undeformed_left_y_Converged = np.copy(mesh_undeformed_left_y_Struct)
mesh_undeformed_left_z_Converged = np.copy(mesh_undeformed_left_z_Struct)

mesh_deformed_left_Converged = np.copy(mesh_deformed_left_Struct)
mesh_deformed_left_x_Converged = np.copy(mesh_deformed_left_x_Struct)
mesh_deformed_left_y_Converged = np.copy(mesh_deformed_left_y_Struct)
mesh_deformed_left_z_Converged = np.copy(mesh_deformed_left_z_Struct)

mesh_delta_left_Converged = np.copy(mesh_delta_left_Struct)
mesh_delta_left_x_Converged = np.copy(mesh_delta_left_x_Struct)
mesh_delta_left_y_Converged = np.copy(mesh_delta_left_y_Struct)
mesh_delta_left_z_Converged = np.copy(mesh_delta_left_z_Struct)

spar_pts_undeformed_left_Converged = np.copy(spar_pts_undeformed_left_Struct)
spar_pts_undeformed_left_x_Converged = np.copy(spar_pts_undeformed_left_x_Struct)
spar_pts_undeformed_left_y_Converged = np.copy(spar_pts_undeformed_left_y_Struct)
spar_pts_undeformed_left_z_Converged = np.copy(spar_pts_undeformed_left_z_Struct)

spar_pts_deformed_left_Converged = np.copy(spar_pts_deformed_left_Struct)
spar_pts_deformed_left_x_Converged = np.copy(spar_pts_deformed_left_x_Struct)
spar_pts_deformed_left_y_Converged = np.copy(spar_pts_deformed_left_y_Struct) 
spar_pts_deformed_left_z_Converged = np.copy(spar_pts_deformed_left_z_Struct)

mesh_point_force_Converged = np.copy(mesh_point_force)

# Print
printSwitch = False
if printSwitch == True:
    print("================================= Results ===================================\n") 
    print("Converged aerodynamics loads:\n", loads_Converged,"\n") 
    print("Converged deformed mesh:\n", mesh_deformed_left_Converged,"\n")  
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
    xpw.PlotWingAndLoads(mesh_initial_left, mesh_undeformed_left_Converged, mesh_deformed_left_Converged,
    spar_pts_undeformed_left_Converged, spar_pts_deformed_left_Converged,
    loads_Converged)

produceThreeJS = True
if produceThreeJS == True:
    fileName = "ThreeJsPlot\models\AeroStruct\AeroStructDesignPoint_AoA_" +  f"{aif.alpha:.1f}" + "_Mach_" + f"{aif.Mach_number:.1f}" + ".json"
    xpt.ThreeJSFile(fileName, mesh_deformed_left_Converged, mesh_point_force_Converged, mesh_delta_left_Converged)

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
    mesh_undeformed_left_Converged, mesh_undeformed_left_x_Converged, mesh_undeformed_left_y_Converged, mesh_undeformed_left_z_Converged,
    spar_pts_undeformed_left_Converged, spar_pts_undeformed_left_x_Converged, spar_pts_undeformed_left_y_Converged, spar_pts_undeformed_left_z_Converged,
    mesh_delta_left_Converged, mesh_delta_left_x_Converged, mesh_delta_left_y_Converged, mesh_delta_left_z_Converged,
    mesh_deformed_left_Converged, mesh_deformed_left_x_Converged, mesh_deformed_left_y_Converged, mesh_deformed_left_z_Converged,
    spar_pts_deformed_left_Converged, spar_pts_deformed_left_x_Converged, spar_pts_deformed_left_y_Converged, spar_pts_deformed_left_z_Converged,
    wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
    loads_Converged, loads_Fx_Converged, loads_Fy_Converged, loads_Fz_Converged, loads_Mx_Converged, loads_My_Converged, loads_Mz_Converged)



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
            
            print("Mach = ", aif.Mach_number, ", alpha = ", aif.alpha)
            # Run analysis
            # iteration to reduce the gap
            
            # reset error
            error = 1
            # reset Counter
            iterCount = 0

            mesh_delta_left_x_Struct_dragPolar = np.copy(mesh_delta_left_x_Converged)
            mesh_delta_left_y_Struct_dragPolar = np.copy(mesh_delta_left_y_Converged)
            mesh_delta_left_z_Struct_dragPolar = np.copy(mesh_delta_left_z_Converged)

            while (error > tol and iterCount < iterMax):

                # pass deformation from Struct to Aero
                mesh_delta_left_x_Aero_dragPolar = np.copy(mesh_delta_left_x_Struct_dragPolar)
                mesh_delta_left_y_Aero_dragPolar = np.copy(mesh_delta_left_y_Struct_dragPolar)
                mesh_delta_left_z_Aero_dragPolar = np.copy(mesh_delta_left_z_Struct_dragPolar)

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
                mesh_delta_left_Aero_dragPolar, # mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
                mesh_deformed_left_Aero_dragPolar, mesh_deformed_left_x_Aero_dragPolar, mesh_deformed_left_y_Aero_dragPolar, mesh_deformed_left_z_Aero_dragPolar,
                spar_pts_deformed_left_Aero_dragPolar, spar_pts_deformed_left_x_Aero_dragPolar, spar_pts_deformed_left_y_Aero_dragPolar, spar_pts_deformed_left_z_Aero_dragPolar,
                wing_area_asref_dragPolar, L_Wing_dragPolar, D_Wing_dragPolar, LoD_Wing_dragPolar,
                CL_Wing_total_dragPolar, CD_Wing_i_dragPolar, CD_Wing_v_dragPolar, CD_Wing_w_dragPolar, CD_Wing_total_dragPolar, 
                CM_Wing_dragPolar, CM_Wing_roll_dragPolar, CM_Wing_pitch_dragPolar, CM_Wing_yaw_dragPolar,
                sec_forces_dragPolar, sec_forces_Fx_dragPolar, sec_forces_Fy_dragPolar, sec_forces_Fz_dragPolar, mesh_point_force_dragPolar,
                loads_Aero_dragPolar, loads_Fx_Aero_dragPolar, loads_Fy_Aero_dragPolar, loads_Fz_Aero_dragPolar, loads_Mx_Aero_dragPolar, loads_My_Aero_dragPolar, loads_Mz_Aero_dragPolar
                ) = aiw.AerodynamicsFunction_wp_10_DM(aif.wing_area_wpref, aif.aspect_ratio_wpref, aif.kink_location_ratio, aif.body_side_ratio,   # Planform definition  
                aif.taper_ratio_trap, aif.root_chord_extension_ratio, aif.trap_quarter_sweep,   # Planform definition
                aif.dihedral, aif.twist_cp_1, aif.twist_cp_2, aif.twist_cp_3, aif.t_over_c_cp_1, aif.t_over_c_cp_2, aif.t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
                aif.CL0, aif.CD0, aif.k_lam, aif.c_max_t, aif.fem_origin,   # Assumed variables
                aif.velocity, aif.alpha, aif.Mach_number, aif.Re, aif.rho,  # Flight conditions
                aif.cg_location_x, aif.cg_location_y, aif.cg_location_z,    # cg locations (if known)
                mesh_delta_left_x_Aero_dragPolar, mesh_delta_left_y_Aero_dragPolar, mesh_delta_left_z_Aero_dragPolar)

                # Pass loads from Aero to Struct
                loads_Struct_dragPolar = np.copy(loads_Aero_dragPolar)

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
                mesh_delta_left_Struct_dragPolar, mesh_delta_left_x_Struct_dragPolar, mesh_delta_left_y_Struct_dragPolar, mesh_delta_left_z_Struct_dragPolar,
                mesh_deformed_left_Struct_dragPolar, mesh_deformed_left_x_Struct_dragPolar, mesh_deformed_left_y_Struct_dragPolar, mesh_deformed_left_z_Struct_dragPolar,
                spar_pts_deformed_left_Struct_dragPolar, spar_pts_deformed_left_x_Struct_dragPolar, spar_pts_deformed_left_y_Struct_dragPolar, spar_pts_deformed_left_z_Struct_dragPolar,
                sparThickness, sparRadius, thickness_intersects, 
                element_mass, wingStructuralMass, vonmisesStress_dragPolar, failure_dragPolar
                ) = siw.StructFunction_wp_10(
                sif.wing_area_wpref, sif.aspect_ratio_wpref, sif.kink_location_ratio, sif.body_side_ratio, sif.taper_ratio_trap, sif.root_chord_extension_ratio, sif.trap_quarter_sweep,   # Planform definition
                sif.dihedral, sif.twist_cp_1, sif.twist_cp_2, sif.twist_cp_3, sif.t_over_c_cp_1, sif.t_over_c_cp_2, sif.t_over_c_cp_3, sif.c_max_t,  # 3D geometry: 1, 2, 3 are from tip to root
                sif.thickness_cp_1, sif.thickness_cp_2, sif.thickness_cp_3, #thickness control point
                sif.E, sif.G, sif.yieldStress, sif.mrho, sif.fem_origin, sif.wing_weight_ratio, # Assumed variables
                loads_Struct_dragPolar # sif.loads_Fx, sif.loads_Fy, sif.loads_Fz, sif.loads_Mx, sif.loads_My, sif.loads_Mz # Loads
                )

                # current gap between the deformed mesh from aero and structure
                error = np.square(np.subtract(mesh_deformed_left_Struct_dragPolar, mesh_deformed_left_Aero_dragPolar)).mean()

                # Count the number of iterations
                iterCount += 1

                print(iterCount,": ",error)
                
            # While loop ends here

            # Converged converged values
            loads_dragPolar_Converged = np.copy(loads_Aero_dragPolar)
            loads_Fx_dragPolar_Converged = np.copy(loads_Fx_Aero_dragPolar)
            loads_Fy_dragPolar_Converged = np.copy(loads_Fy_Aero_dragPolar)
            loads_Fz_dragPolar_Converged = np.copy(loads_Fz_Aero_dragPolar)
            loads_Mx_dragPolar_Converged = np.copy(loads_Mx_Aero_dragPolar)
            loads_My_dragPolar_Converged = np.copy(loads_My_Aero_dragPolar)
            loads_Mz_dragPolar_Converged = np.copy(loads_Mz_Aero_dragPolar)

            # mesh_undeformed_left_Converged = np.copy(mesh_undeformed_left_Struct)
            # mesh_undeformed_left_x_Converged = np.copy(mesh_undeformed_left_x_Struct)
            # mesh_undeformed_left_y_Converged = np.copy(mesh_undeformed_left_y_Struct)
            # mesh_undeformed_left_z_Converged = np.copy(mesh_undeformed_left_z_Struct)

            mesh_deformed_left_dragPolar_Converged = np.copy(mesh_deformed_left_Struct_dragPolar)
            mesh_deformed_left_x_dragPolar_Converged = np.copy(mesh_deformed_left_x_Struct_dragPolar)
            mesh_deformed_left_y_dragPolar_Converged = np.copy(mesh_deformed_left_y_Struct_dragPolar)
            mesh_deformed_left_z_dragPolar_Converged = np.copy(mesh_deformed_left_z_Struct_dragPolar)

            mesh_delta_left_dragPolar_Converged = np.copy(mesh_delta_left_Struct_dragPolar)
            mesh_delta_left_x_dragPolar_Converged = np.copy(mesh_delta_left_x_Struct_dragPolar)
            mesh_delta_left_y_dragPolar_Converged = np.copy(mesh_delta_left_y_Struct_dragPolar)
            mesh_delta_left_z_dragPolar_Converged = np.copy(mesh_delta_left_z_Struct_dragPolar)

            spar_pts_undeformed_left_Converged = np.copy(spar_pts_undeformed_left_Struct)
            spar_pts_undeformed_left_x_Converged = np.copy(spar_pts_undeformed_left_x_Struct)
            spar_pts_undeformed_left_y_Converged = np.copy(spar_pts_undeformed_left_y_Struct)
            spar_pts_undeformed_left_z_Converged = np.copy(spar_pts_undeformed_left_z_Struct)

            spar_pts_deformed_left_dragPolar_Converged = np.copy(spar_pts_deformed_left_Struct_dragPolar)
            spar_pts_deformed_left_x_dragPolar_Converged = np.copy(spar_pts_deformed_left_x_Struct_dragPolar)
            spar_pts_deformed_left_y_dragPolar_Converged = np.copy(spar_pts_deformed_left_y_Struct_dragPolar) 
            spar_pts_deformed_left_z_dragPolar_Converged = np.copy(spar_pts_deformed_left_z_Struct_dragPolar)

            CL_Wing_total_dragPolar_Converged = CL_Wing_total_dragPolar
            CD_Wing_total_dragPolar_Converged = CD_Wing_total_dragPolar
            CD_Wing_i_dragPolar_Converged = CD_Wing_i_dragPolar
            mesh_point_force_dragPolar_Converged = np.copy(mesh_point_force_dragPolar)

            # Record CL, CD
            CL_dragPolar_Raw[i,j] = CL_Wing_total_dragPolar
            CD_total_dragPolar_Raw[i,j] = CD_Wing_total_dragPolar
            CD_i_dragPolar_raw[i,j] = CD_Wing_i_dragPolar # induced drag

            produceThreeJS = True
            if produceThreeJS == True:
                fileName = "ThreeJsPlot\models\AeroStruct\DragPolar\AeroStructDragPolar_AoA_" +  f"{alpha_dragPolar[j]:.1f}" + "_Mach_" + f"{Mach_dragPolar[i]:.1f}" + ".json"
                xpt.ThreeJSFile(fileName, mesh_deformed_left_dragPolar_Converged, mesh_point_force_dragPolar_Converged, mesh_delta_left_dragPolar_Converged)

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