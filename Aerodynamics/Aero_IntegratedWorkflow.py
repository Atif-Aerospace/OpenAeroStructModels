import numpy as np

import Common.Common_PlanformRelationship as cpr
import Common.Common_InitialMesh as cim
import Aerodynamics.Aero_UnderformedMesh as aum
import Aerodynamics.Aero_DerformedMesh as adm
import Aerodynamics.Aero_AnalysisComponent as aac

# Define the wing with gross wing area Option 1 and leading edge sweep option 1
# Deformation defined with control points
def AerodynamicsFunction_geo_11_CP(wing_area_geo, aspect_ratio_geo, kink_location_ratio, taper_ratio_t_r, taper_ratio_k_r, inboard_LE_sweep, outboard_LE_sweep,   # Planform definition
    dihedral, twist_cp_1, twist_cp_2, twist_cp_3, t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    CL0, CD0, k_lam, c_max_t, fem_origin,   # Assumed variables
    velocity, alpha, Mach_number, Re, rho,  # Flight conditions
    cg_location_x, cg_location_y, cg_location_z,    # cg locations (if known)
    delta_x_LE_1, delta_x_LE_2, delta_x_LE_3,   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
    delta_y_LE_1, delta_y_LE_2, delta_y_LE_3,
    delta_z_LE_1, delta_z_LE_2, delta_z_LE_3,
    delta_x_TE_1, delta_x_TE_2, delta_x_TE_3,   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
    delta_y_TE_1, delta_y_TE_2, delta_y_TE_3,
    delta_z_TE_1, delta_z_TE_2, delta_z_TE_3):

    # Conversion of related planform geometry parameters
    
    # Sizing
    (entire_span, half_span, kink_location, root_chord, kink_chord, tip_chord, 
    wingarea_dependency, aspectratio_dependency, entirespan_dependency, 
    halfspan_dependency, kinklocation_dependency, kinklocation_ratio_dependency,
    taperratioTr_dependency, taperratioKr_dependency, 
    rootchord_dependency, kinkchord_dependency, tipchord_dependency
    ) = cpr.SizeRelation_1(wing_area_geo, aspect_ratio_geo, kink_location_ratio, taper_ratio_t_r, taper_ratio_k_r)
    
    # Sweep angle
    (inboard_quarter_sweep, outboard_quarter_sweep, 
    inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency
    ) = cpr.SweepRelation_1(inboard_LE_sweep, outboard_LE_sweep, half_span, kink_location, root_chord, kink_chord, tip_chord)
    
    # Define initial mesh, No twist and dihedral, No structural deformation
    # nx = 5  # number of chordwise nodal points (should be odd)
    # ny_outboard = 12  # number of spanwise nodal points for the outboard segment
    # ny_inboard = 5  # number of spanwise nodal points for the inboard segment
    
    ny = 26
    ny_inboard = round(kink_location_ratio * ny)
    ny_outboard = ny + 1 - ny_inboard
    # nx = round(root_chord / half_span * ny)
    nx = 7

    (ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z
    ) = cim.InitialMesh(half_span, kink_location, root_chord, kink_chord, tip_chord, inboard_LE_sweep, outboard_LE_sweep, nx, ny_outboard, ny_inboard)

    # mesh after twist and dihedral, no structural deformation yet
    twist_cp = np.array([twist_cp_1, twist_cp_2, twist_cp_3])
    t_over_c_cp = np.array([t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3])

    (t_over_c_actual, twist_actual, wingInfo,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z
    ) = aum.UndeformedMesh(mesh_initial_right, dihedral, twist_cp, t_over_c_cp, CL0, CD0, k_lam, c_max_t, fem_origin)

    # mesh_delta: displacement at each mesh point due to structural deformation, 
    # in this aerodynamic analysis, this is either from the results of structural analysis or manually added by using the spline line 
    # the latter can be used for model reversal
    (mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z
    ) = adm.DeformedMesh_CP(mesh_undeformed_left, fem_origin,
                        delta_x_LE_1, delta_x_LE_2, delta_x_LE_3,   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
                        delta_y_LE_1, delta_y_LE_2, delta_y_LE_3,
                        delta_z_LE_1, delta_z_LE_2, delta_z_LE_3,
                        delta_x_TE_1, delta_x_TE_2, delta_x_TE_3,   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
                        delta_y_TE_1, delta_y_TE_2, delta_y_TE_3,
                        delta_z_TE_1, delta_z_TE_2, delta_z_TE_3)

    # run the actual aerodynamics analysis using OpenAeroStruct components
    cg_location = np.array([cg_location_x, cg_location_y, cg_location_z])
    (wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, 
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
    ) = aac.AerodynamicsAnalysis(wingInfo, velocity, alpha, Mach_number, Re, rho, cg_location, mesh_deformed_left, t_over_c_actual, mesh_delta_left)

    return (entire_span, half_span, kink_location, root_chord, kink_chord, tip_chord, inboard_quarter_sweep, outboard_quarter_sweep, 
    wingarea_dependency, aspectratio_dependency, entirespan_dependency, 
    halfspan_dependency, kinklocation_dependency, kinklocation_ratio_dependency,
    taperratioTr_dependency, taperratioKr_dependency, 
    rootchord_dependency, kinkchord_dependency, tipchord_dependency,
    inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency,
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
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz)

# Define the wing with wimpress reference area Option 1 and quarter chord sweep option 0
# Deformation defined with control points
def AerodynamicsFunction_wp_10_CP(wing_area_wpref, aspect_ratio_wpref, kink_location_ratio, body_side_ratio, taper_ratio_trap, root_chord_extension_ratio, trap_quarter_sweep,   # Planform definition
    dihedral, twist_cp_1, twist_cp_2, twist_cp_3, t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    CL0, CD0, k_lam, c_max_t, fem_origin,   # Assumed variables
    velocity, alpha, Mach_number, Re, rho,  # Flight conditions
    cg_location_x, cg_location_y, cg_location_z,    # cg locations (if known)
    delta_x_LE_1, delta_x_LE_2, delta_x_LE_3,   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
    delta_y_LE_1, delta_y_LE_2, delta_y_LE_3,
    delta_z_LE_1, delta_z_LE_2, delta_z_LE_3,
    delta_x_TE_1, delta_x_TE_2, delta_x_TE_3,   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
    delta_y_TE_1, delta_y_TE_2, delta_y_TE_3,
    delta_z_TE_1, delta_z_TE_2, delta_z_TE_3):

    # Conversion of related planform geometry parameters
    
    # Sizing
    (wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap, entire_span, half_span, kink_location, 
    root_chord_trap, root_chord, kink_chord, tip_chord, taper_ratio_t_r, taper_ratio_k_r,
    wingarea_wpref_dependency, aspectratio_wpref_dependency, 
    wingarea_geo_dependency, aspectratio_geo_dependency, 
    wingarea_trap_dependency, aspectratio_trap_dependency,
    kinklocation_ratio_dependency, bodysideratio_dependency, rootchord_extratio_dependency, taperratioTrap_dependency, taperratioTr_dependency, taperratioKr_dependency,
    entirespan_dependency, halfspan_dependency, kinklocation_dependency, rootchordTrap_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency
    ) = cpr.SizeRelation_Wimpress_1(wing_area_wpref, aspect_ratio_wpref, kink_location_ratio, body_side_ratio, taper_ratio_trap, root_chord_extension_ratio)
    
    # Sweep angle
    (inboard_quarter_sweep, outboard_quarter_sweep, inboard_LE_sweep, outboard_LE_sweep, 
    trapquartersweep_dependency, inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency
    ) = cpr.SweepRelation_Wimpress_0(trap_quarter_sweep, kink_location, root_chord_trap, kink_chord, root_chord)
    
    # Define initial mesh, No twist and dihedral, No structural deformation
    # nx = 5  # number of chordwise nodal points (should be odd)
    # ny_outboard = 12  # number of spanwise nodal points for the outboard segment
    # ny_inboard = 5  # number of spanwise nodal points for the inboard segment
    
    ny = 15
    ny_inboard = round(kink_location_ratio * ny)
    # ny_inboard = 3
    ny_outboard = ny + 1 - ny_inboard
    # nx = round(root_chord / half_span * ny)
    nx = 7

    (ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z
    ) = cim.InitialMesh(half_span, kink_location, root_chord, kink_chord, tip_chord, inboard_LE_sweep, outboard_LE_sweep, nx, ny_outboard, ny_inboard)

    # mesh after twist and dihedral, no structural deformation yet
    twist_cp = np.array([twist_cp_1, twist_cp_2, twist_cp_3])
    t_over_c_cp = np.array([t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3])

    (t_over_c_actual, twist_actual, wingInfo,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z
    ) = aum.UndeformedMesh(mesh_initial_right, dihedral, twist_cp, t_over_c_cp, CL0, CD0, k_lam, c_max_t, fem_origin)

    # mesh_delta: displacement at each mesh point due to structural deformation, 
    # in this aerodynamic analysis, this is either from the results of structural analysis or manually added by using the spline line 
    # the latter can be used for model reversal
    (mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z
    ) = adm.DeformedMesh_CP(mesh_undeformed_left, fem_origin,
                        delta_x_LE_1, delta_x_LE_2, delta_x_LE_3,   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
                        delta_y_LE_1, delta_y_LE_2, delta_y_LE_3,
                        delta_z_LE_1, delta_z_LE_2, delta_z_LE_3,
                        delta_x_TE_1, delta_x_TE_2, delta_x_TE_3,   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
                        delta_y_TE_1, delta_y_TE_2, delta_y_TE_3,
                        delta_z_TE_1, delta_z_TE_2, delta_z_TE_3)

    # run the actual aerodynamics analysis using OpenAeroStruct components
    cg_location = np.array([cg_location_x, cg_location_y, cg_location_z])
    (wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total,
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
    ) = aac.AerodynamicsAnalysis(wingInfo, velocity, alpha, Mach_number, Re, rho, cg_location, mesh_deformed_left, t_over_c_actual, mesh_delta_left)

    return (wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap, entire_span, half_span, kink_location,
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
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz)


# Define the wing with wimpress reference area Option 1 and quarter chord sweep option 0
# Define the defromed mesh based on mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z
def AerodynamicsFunction_wp_10_DM(wing_area_wpref, aspect_ratio_wpref, kink_location_ratio, body_side_ratio, taper_ratio_trap, root_chord_extension_ratio, trap_quarter_sweep,   # Planform definition
    dihedral, twist_cp_1, twist_cp_2, twist_cp_3, t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    CL0, CD0, k_lam, c_max_t, fem_origin,   # Assumed variables
    velocity, alpha, Mach_number, Re, rho,  # Flight conditions
    cg_location_x, cg_location_y, cg_location_z,    # cg locations (if known)
    mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z):

    # Conversion of related planform geometry parameters
    
    # Sizing
    (wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap, entire_span, half_span, kink_location, 
    root_chord_trap, root_chord, kink_chord, tip_chord, taper_ratio_t_r, taper_ratio_k_r,
    wingarea_wpref_dependency, aspectratio_wpref_dependency, 
    wingarea_geo_dependency, aspectratio_geo_dependency, 
    wingarea_trap_dependency, aspectratio_trap_dependency,
    kinklocation_ratio_dependency, bodysideratio_dependency, rootchord_extratio_dependency, taperratioTrap_dependency, taperratioTr_dependency, taperratioKr_dependency,
    entirespan_dependency, halfspan_dependency, kinklocation_dependency, rootchordTrap_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency
    ) = cpr.SizeRelation_Wimpress_1(wing_area_wpref, aspect_ratio_wpref, kink_location_ratio, body_side_ratio, taper_ratio_trap, root_chord_extension_ratio)
    
    # Sweep angle
    (inboard_quarter_sweep, outboard_quarter_sweep, inboard_LE_sweep, outboard_LE_sweep, 
    trapquartersweep_dependency, inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency
    ) = cpr.SweepRelation_Wimpress_0(trap_quarter_sweep, kink_location, root_chord_trap, kink_chord, root_chord)
    
    # Define initial mesh, No twist and dihedral, No structural deformation
    # nx = 5  # number of chordwise nodal points (should be odd)
    # ny_outboard = 12  # number of spanwise nodal points for the outboard segment
    # ny_inboard = 5  # number of spanwise nodal points for the inboard segment
    
    ny = 15
    ny_inboard = round(kink_location_ratio * ny)
    # ny_inboard = 3
    ny_outboard = ny + 1 - ny_inboard
    # nx = round(root_chord / half_span * ny)
    nx = 7

    (ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z
    ) = cim.InitialMesh(half_span, kink_location, root_chord, kink_chord, tip_chord, inboard_LE_sweep, outboard_LE_sweep, nx, ny_outboard, ny_inboard)

    # mesh after twist and dihedral, no structural deformation yet
    twist_cp = np.array([twist_cp_1, twist_cp_2, twist_cp_3])
    t_over_c_cp = np.array([t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3])

    (t_over_c_actual, twist_actual, wingInfo,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z
    ) = aum.UndeformedMesh(mesh_initial_right, dihedral, twist_cp, t_over_c_cp, CL0, CD0, k_lam, c_max_t, fem_origin)

    # mesh_delta: displacement at each mesh point due to structural deformation, 
    # in this aerodynamic analysis, this is either from the results of structural analysis or manually added by using the spline line 
    # the latter can be used for model reversal
    (mesh_delta_left, # mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z
    ) = adm.DeformedMesh_DM(mesh_undeformed_left, fem_origin, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z)

    # run the actual aerodynamics analysis using OpenAeroStruct components
    cg_location = np.array([cg_location_x, cg_location_y, cg_location_z])
    (wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total,
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
    ) = aac.AerodynamicsAnalysis(wingInfo, velocity, alpha, Mach_number, Re, rho, cg_location, mesh_deformed_left, t_over_c_actual, mesh_delta_left) # mesh_delta_left is used for three.js plot

    return (wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap, entire_span, half_span, kink_location,
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
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz)