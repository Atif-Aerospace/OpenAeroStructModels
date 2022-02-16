import numpy as np

import Common.Common_PlanformRelationship as cpr
import Common.Common_InitialMesh as cim
import Structures.Struct_UnderformedMesh as sum
import Structures.Struct_AnalysisComponent as sac

def StructFunction_wp_10( 
    wing_area_wpref, aspect_ratio_wpref, kink_location_ratio, body_side_ratio, taper_ratio_trap, root_chord_extension_ratio, trap_quarter_sweep,   # Planform definition
    dihedral, twist_cp_1, twist_cp_2, twist_cp_3, t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3, c_max_t,  # 3D geometry: 1, 2, 3 are from tip to root
    thickness_cp_1, thickness_cp_2, thickness_cp_3, #thickness control point
    E, G, yieldStress, mrho, fem_origin, wing_weight_ratio, # Assumed variables
    loads # loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
    ):

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
    ny_outboard = ny + 1 - ny_inboard
    # nx = round(root_chord / half_span * ny)
    nx = 7

    (ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z
    ) = cim.InitialMesh(half_span, kink_location, root_chord, kink_chord, tip_chord, inboard_LE_sweep, outboard_LE_sweep, nx, ny_outboard, ny_inboard)

    # mesh after twist and dihedral, no structural deformation yet
    twist_cp = np.array([twist_cp_1, twist_cp_2, twist_cp_3])
    t_over_c_cp = np.array([t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3])
    thickness_cp = np.array([thickness_cp_1, thickness_cp_2, thickness_cp_3])

    (t_over_c_actual, twist_actual, wingInfo,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z
    ) = sum.UndeformedMesh(mesh_initial_right, dihedral, twist_cp, thickness_cp, t_over_c_cp, c_max_t, E, G, yieldStress, mrho, fem_origin, wing_weight_ratio)


    # ny = wingInfo["mesh"].shape[1]
    # loads = np.ones((ny, 6)) * 1e5
    # loads[:,0] = loads_Fx
    # loads[:,1] = loads_Fy
    # loads[:,2] = loads_Fz
    # loads[:,3] = loads_Mx
    # loads[:,4] = loads_My
    # loads[:,5] = loads_Mz


    # run the actual structural analysis using OpenAeroStruct components
    (mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
        mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z, 
        spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
        sparThickness, sparRadius, thickness_intersects, 
        element_mass, wingStructuralMass, vonmisesStress, failure
    ) = sac.StructualAnalysis(mesh_undeformed_left, wingInfo, loads)




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
    t_over_c_actual, twist_actual, wingInfo, loads,
    mesh_undeformed_left, mesh_undeformed_left_x, mesh_undeformed_left_y, mesh_undeformed_left_z,
    spar_pts_undeformed_left, spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z,
    mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
    mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
    spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z,
    sparThickness, sparRadius, thickness_intersects, 
    element_mass, wingStructuralMass, vonmisesStress, failure)