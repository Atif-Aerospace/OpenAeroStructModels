import numpy as np

def WriteToFiles(fileName, wing_area_wpref, aspect_ratio_wpref, wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap,  # Planform definition 
    kink_location_ratio, body_side_ratio, taper_ratio_t_r, taper_ratio_k_r, taper_ratio_trap, root_chord_extension_ratio,    # Planform definition
    entire_span, half_span, kink_location, root_chord_trap, root_chord, kink_chord, tip_chord, 
    inboard_quarter_sweep, outboard_quarter_sweep, inboard_LE_sweep, outboard_LE_sweep, trap_quarter_sweep,   # Planform definition
    dihedral, twist_cp_1, twist_cp_2, twist_cp_3, t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
    CL0, CD0, k_lam, c_max_t, fem_origin,   # Assumed variables
    velocity, alpha, Mach_number, Re, rho,  # Flight conditions
    cg_location_x, cg_location_y, cg_location_z,    # cg locations (if known)
    delta_x_LE_1, delta_x_LE_2, delta_x_LE_3,   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
    delta_y_LE_1, delta_y_LE_2, delta_y_LE_3,
    delta_z_LE_1, delta_z_LE_2, delta_z_LE_3,
    delta_x_TE_1, delta_x_TE_2, delta_x_TE_3,   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
    delta_y_TE_1, delta_y_TE_2, delta_y_TE_3,
    delta_z_TE_1, delta_z_TE_2, delta_z_TE_3,
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
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz):


    f = open(fileName, "w")

    f.write("=====================================================================================\n")
    f.write("||                               Wing Definition Summary                           ||\n")
    f.write("=====================================================================================\n")

    f.write("\n================================ Planform Specification =============================\n")
    f.write("\nwing_area_trap               = ")
    f.write("{:.6f}".format(wing_area_trap))
    f.write("     " + wingarea_trap_dependency)
    f.write("       (Trapezoid part of the wing in Wimpress defintion)")

    f.write("\nwing_area_wpref              = ")
    f.write("{:.6f}".format(wing_area_wpref))
    f.write("     " + wingarea_wpref_dependency)
    f.write("       (Reference area based on Wimpress defintion)")

    f.write("\nwing_area_geo                = ")
    f.write("{:.6f}".format(wing_area_geo))
    f.write("     " + wingarea_geo_dependency)
    f.write("       (Gross area based on the entire geometry)")

    f.write("\nwing_area_asref              = ")
    f.write("{:.6f}".format(wing_area_asref))
    f.write("     (Dependent)       (Computed by OpenAeroStruct VLMGeometry component)")

    f.write("\naspect_ratio_trap            = ")
    f.write("{:.6f}".format(aspect_ratio_trap))
    f.write("       " + aspectratio_trap_dependency)
    f.write("       (Aspect ratio of the trapezoid part in Wimpress defintion)")

    f.write("\naspect_ratio_wpref           = ")
    f.write("{:.6f}".format(aspect_ratio_wpref))
    f.write("       " + aspectratio_wpref_dependency)
    f.write("       (Aspect ratio based on Wimpress reference area)")

    f.write("\naspect_ratio_geo             = ")
    f.write("{:.6f}".format(aspect_ratio_geo))
    f.write("       " + aspectratio_geo_dependency)

    f.write("\nentire_span                  = ")
    f.write("{:.6f}".format(entire_span))
    f.write("      " + entirespan_dependency)

    f.write("\nhalf_span                    = ")
    f.write("{:.6f}".format(half_span))
    f.write("      " + halfspan_dependency)

    f.write("\nbody_side_ratio              = ")
    f.write("{:.6f}".format(body_side_ratio))
    f.write("       " + bodysideratio_dependency)

    f.write("\nkink_location_ratio          = ")
    f.write("{:.6f}".format(kink_location_ratio))
    f.write("       " + kinklocation_ratio_dependency)

    f.write("\nkink_location                = ")
    f.write("{:.6f}".format(kink_location))
    f.write("      " + kinklocation_dependency)

    f.write("\ntaper_ratio_trap             = ")
    f.write("{:.6f}".format(taper_ratio_trap))
    f.write("       " + taperratioTrap_dependency)
    f.write("       (Taper ratio of the trapezoid part in Wimpress defintion)")

    f.write("\ntaper_ratio_t_r              = ")
    f.write("{:.6f}".format(taper_ratio_t_r))
    f.write("       " + taperratioTr_dependency)

    f.write("\ntaper_ratio_k_r              = ")
    f.write("{:.6f}".format(taper_ratio_k_r))
    f.write("       " + taperratioKr_dependency)

    f.write("\nroot_chord_extension_ratio   = ")
    f.write("{:.6f}".format(root_chord_extension_ratio))
    f.write("       " + rootchord_extratio_dependency)
    f.write("       (Actual root chord / Root chord of the trapezoid part)")

    f.write("\nroot_chord_trap              = ")
    f.write("{:.6f}".format(root_chord_trap))
    f.write("       " + rootchordTrap_dependency)
    f.write("       (Root chord of the trapezoid part in Wimpress defintion)")

    f.write("\nroot_chord                   = ")
    f.write("{:.6f}".format(root_chord))
    f.write("      " + rootchord_dependency)

    f.write("\nkink_chord                   = ")
    f.write("{:.6f}".format(kink_chord))
    f.write("       " + kinkchord_dependency)

    f.write("\ntip_chord                    = ")
    f.write("{:.6f}".format(tip_chord))
    f.write("       " + tipchord_dependency)

    f.write("\ntrap_quarter_sweep           = ")
    f.write("{:.6f}".format(trap_quarter_sweep))
    f.write("      " + trapquartersweep_dependency)

    f.write("\ninboard_quarter_sweep        = ")
    f.write("{:.6f}".format(inboard_quarter_sweep))
    f.write("      " + inquartersweep_dependency)

    f.write("\noutboard_quarter_sweep       = ")
    f.write("{:.6f}".format(outboard_quarter_sweep))
    f.write("      " + outquartersweep_dependency)

    f.write("\ninboard_LE_sweep             = ")
    f.write("{:.6f}".format(inboard_LE_sweep))
    f.write("      " + inLEsweep_dependency)

    f.write("\noutboard_LE_sweep            = ")
    f.write("{:.6f}".format(outboard_LE_sweep))
    f.write("      " + outLEsweep_dependency)

    f.write("\n\n============================ 3D allocation of the sections ==========================\n")

    twist_cp = np.array([twist_cp_1, twist_cp_2, twist_cp_3])
    f.write("\ntwist_cp (from tip to root):\n")
    twist_cp_str = np.array2string(twist_cp, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(twist_cp_str + "\n")

    f.write("\nActual_twist [deg] (from tip to root):\n")
    twist_actual_str = np.array2string(twist_actual, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(twist_actual_str + "\n")

    t_over_c_cp = np.array([t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3])
    f.write("\nt_over_c_cp (from tip to root):\n")
    t_over_c_cp_str = np.array2string(t_over_c_cp, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(t_over_c_cp_str + "\n")

    f.write("\nActual_t_over_c (from tip to root):\n")
    t_over_c_str = np.array2string(t_over_c_actual, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(t_over_c_str + "\n")

    f.write("\n\n=================================== Initial Mesh =====================================\n")
    f.write("this is the mesh (from tip to root) from the wing planform definition (no twist and dihedral)\n")
    f.write("\nmesh_initial_left_x:\n")
    mesh_initial_left_x_str = np.array2string(mesh_initial_left_x, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_initial_left_x_str + "\n")

    f.write("\nmesh_initial_left_y:\n")
    mesh_initial_left_y_str = np.array2string(mesh_initial_left_y, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_initial_left_y_str + "\n")

    f.write("\nmesh_initial_left_z:\n")
    mesh_initial_left_z_str = np.array2string(mesh_initial_left_z, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_initial_left_z_str + "\n")

    f.write("\n\n============================== Undeformed Mesh =====================================\n")
    f.write("Mesh (from tip to root) after twist and dihedral?, no structural deformation\n")
    f.write("\nmesh_undeformed_left_x:\n")
    mesh_undeformed_left_x_str = np.array2string(mesh_undeformed_left_x, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_undeformed_left_x_str + "\n")

    f.write("\nmesh_undeformed_left_y:\n")
    mesh_undeformed_left_y_str = np.array2string(mesh_undeformed_left_y, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_undeformed_left_y_str + "\n")

    f.write("\nmesh_undeformed_left_z:\n")
    mesh_undeformed_left_z_str = np.array2string(mesh_undeformed_left_z, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_undeformed_left_z_str + "\n")


    f.write("\n\n========================= Displacement at Each Mesh Point ============================\n")

    f.write("\ndisplacement pre-defined at three points on the leading edge, from TIP to ROOT\n")

    delta_x_LE = np.array([delta_x_LE_1, delta_x_LE_2, delta_x_LE_3])
    f.write("\ndelta_x_LE:\n")
    delta_x_LE_str = np.array2string(delta_x_LE, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(delta_x_LE_str + "\n")

    delta_y_LE = np.array([delta_y_LE_1, delta_y_LE_2, delta_y_LE_3])
    f.write("\ndelta_y_LE:\n")
    delta_y_LE_str = np.array2string(delta_y_LE, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(delta_y_LE_str + "\n")

    delta_z_LE = np.array([delta_z_LE_1, delta_z_LE_2, delta_z_LE_3])
    f.write("\ndelta_z_LE:\n")
    delta_z_LE_str = np.array2string(delta_z_LE, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(delta_z_LE_str + "\n")

    f.write("\ndisplacement pre-defined at three points on the trailing edge, from TIP to ROOT\n")

    delta_x_TE = np.array([delta_x_TE_1, delta_x_TE_2, delta_x_TE_3])
    f.write("\ndelta_x_TE:\n")
    delta_x_TE_str = np.array2string(delta_x_TE, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(delta_x_TE_str + "\n")

    delta_y_TE = np.array([delta_y_TE_1, delta_y_TE_2, delta_y_TE_3])
    f.write("\ndelta_y_TE:\n")
    delta_y_TE_str = np.array2string(delta_y_TE, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(delta_y_TE_str + "\n")

    delta_z_TE = np.array([delta_z_TE_1, delta_z_TE_2, delta_z_TE_3])
    f.write("\ndelta_z_TE:\n")
    delta_z_TE_str = np.array2string(delta_z_TE, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(delta_z_TE_str + "\n")

    f.write("displacement at each mesh point (from tip to root) due to structural deformation, in this aerodynamic analysis,\nthis is either from the results of structural analysis or manually added by using the spline line\nthe latter can be used for model reversal\n")

    f.write("\nmesh_delta_x:\n")
    mesh_delta_left_x_str = np.array2string(mesh_delta_left_x, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_delta_left_x_str + "\n")

    f.write("\nmesh_delta_y:\n")
    mesh_delta_left_y_str = np.array2string(mesh_delta_left_y, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_delta_left_y_str + "\n")

    f.write("\nmesh_delta_z:\n")
    mesh_delta_left_z_str = np.array2string(mesh_delta_left_z, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_delta_left_z_str + "\n")


    f.write("\n\n======================= Actual Mesh Used in Aerodynamics Analysis ===================\n")
    f.write("\nmesh (from tip to root) after the structural deformation, this mesh will be used in the actual aerodynamics analysis\n")
    f.write("should be the same as input_var.def_mesh, aeropoint_group.wing.def_mesh, aeropoint_group.aero_states.wing_def_mesh, and loadtransfer_group.def_mesh\n")
    f.write("\nmesh_deformed_x:\n")
    mesh_deformed_left_x_str = np.array2string(mesh_deformed_left_x, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_deformed_left_x_str + "\n")

    f.write("\nmesh_deformed_y:\n")
    mesh_deformed_left_y_str = np.array2string(mesh_deformed_left_y, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_deformed_left_y_str + "\n")

    f.write("\nmesh_deformed_z:\n")
    mesh_deformed_left_z_str = np.array2string(mesh_deformed_left_z, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(mesh_deformed_left_z_str + "\n")

    f.write("\n\n================================== Assumed Variables ================================\n")

    f.write("\nCL0                      = ")
    f.write("{:.6f}".format(CL0))

    f.write("\nCD0                      = ")
    f.write("{:.6f}".format(CD0))

    f.write("\nk_lam                    = ")
    f.write("{:.6f}".format(k_lam))

    f.write("\nc_max_t                  = ")
    f.write("{:.6f}".format(c_max_t))

    f.write("\nfem_origin               = ")
    f.write("{:.6f}".format(fem_origin))


    f.write("\n\n=====================================================================================\n")
    f.write("||                            Aerodynamic Analysis Results                         ||\n")
    f.write("=====================================================================================\n")

    f.write("\n=============================== Lift and Drag ======================================\n")

    f.write("\nL_Wing                   = ")
    f.write("{:.6f}".format(L_Wing))

    f.write("\nD_Wing                   = ")
    f.write("{:.6f}".format(D_Wing))

    f.write("\n\n================================ Coefficients ======================================\n")

    f.write("\nCL_Wing_total            = ")
    f.write("{:.6f}".format(CL_Wing_total))

    f.write("\nCD_Wing_induced          = ")
    f.write("{:.6f}".format(CD_Wing_i))

    f.write("\nCD_Wing_viscous          = ")
    f.write("{:.6f}".format(CD_Wing_v))

    f.write("\nCD_Wing_wave             = ")
    f.write("{:.6f}".format(CD_Wing_w))

    f.write("\nCD_Wing_total            = ")
    f.write("{:.6f}".format(CD_Wing_total))

    f.write("\nCM_Wing_roll             = ")
    f.write("{:.6f}".format(CM_Wing_roll))

    f.write("\nCM_Wing_pitch            = ")
    f.write("{:.6f}".format(CM_Wing_pitch))

    f.write("\nCM_Wing_yaw              = ")
    f.write("{:.6f}".format(CM_Wing_yaw))

    f.write("\nLift_over_Drag           = ")
    f.write("{:.6f}".format(LoD_Wing))

    f.write("\n\n======================================= Loads ======================================\n")

    f.write("\nloads_Fx [N]:\n")
    loads_Fx_str = np.array2string(loads_Fx, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(loads_Fx_str + "\n")

    f.write("\nloads_Fy [N]:\n")
    loads_Fy_str = np.array2string(loads_Fy, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(loads_Fy_str + "\n")

    f.write("\nloads_Fz [N]:\n")
    loads_Fz_str = np.array2string(loads_Fz, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(loads_Fz_str + "\n")

    f.write("\nloads_Mx:\n")
    loads_Mx_str = np.array2string(loads_Mx, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(loads_Mx_str + "\n")

    f.write("\nloads_My:\n")
    loads_My_str = np.array2string(loads_My, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(loads_My_str + "\n")

    f.write("\nloads_Mz:\n")
    loads_Mz_str = np.array2string(loads_Mz, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    f.write(loads_Mz_str + "\n")

    f.close()

