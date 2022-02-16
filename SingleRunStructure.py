import Structures.Struct_IntegratedWorkflow as siw
import Structures.Struct_InputFileWimpress as sif

import Auxiliary.Aux_PlotWingAndLoads as xpw


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
    ) = siw.StructFunction_wp_10(sif.wing_area_wpref, sif.aspect_ratio_wpref, sif.kink_location_ratio, sif.body_side_ratio, sif.taper_ratio_trap, sif.root_chord_extension_ratio, sif.trap_quarter_sweep,   # Planform definition
    sif.dihedral, sif.twist_cp_1, sif.twist_cp_2, sif.twist_cp_3, sif.t_over_c_cp_1, sif.t_over_c_cp_2, sif.t_over_c_cp_3, sif.c_max_t,  # 3D geometry: 1, 2, 3 are from tip to root
    sif.thickness_cp_1, sif.thickness_cp_2, sif.thickness_cp_3, #thickness control point
    sif.E, sif.G, sif.yieldStress, sif.mrho, sif.fem_origin, sif.wing_weight_ratio, # related to structure
    sif.loads # loads have to be assumed: sif.loads_Fx, sif.loads_Fy, sif.loads_Fz, sif.loads_Mx, sif.loads_My, sif.loads_Mz # Loads
    )
    
print("done")


# Plot
plotSwitch = True
if plotSwitch == True:
    xpw.PlotWingAndLoads(mesh_initial_left, mesh_undeformed_left, mesh_deformed_left,
    spar_pts_undeformed_left, spar_pts_deformed_left,
    loads)