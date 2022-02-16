"""this input file is used for validation, the values are changed to match the example from 
http://mdolab.engin.umich.edu/OpenAeroStruct/aero_walkthrough.html

option = 0: results from my code is very close to OpenAeroStruct(OAS), the difference is caused by initial mesh
option = 1: Two cases are identical

"with_wave": False: Tested wo cases are identical
"with_wave": True: Tested wo cases are identical

"""






###########################################################################
#######                       Run the my code                       #######
###########################################################################

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
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz, mesh_point_force,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz
    ) = aiw.AerodynamicsFunction_wp_10_CP(aif.wing_area_wpref, aif.aspect_ratio_wpref, aif.kink_location_ratio, aif.body_side_ratio,   # Planform definition  
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










###########################################################################
#######                 Run the original OpenAeroAtruct code        #######
###########################################################################
import numpy as np
import openmdao.api as om
from openmdao.api import n2
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint

import matplotlib.pyplot as plt


# Create a dictionary to store options about the mesh

# Original setting
# mesh_dict = {"num_y": 7, "num_x": 2, "wing_type": "CRM", "symmetry": True, "num_twist_cp": 5}

# new setting to match my functions
mesh_dict = {"num_y": 13, "num_x": 3, "wing_type": "CRM", "symmetry": True, "num_twist_cp": 3}

option = 1
if option == 0:
    # Generate the aerodynamic mesh based on the previous dictionary
    mesh, twist_cp = generate_mesh(mesh_dict)
else:
    mesh = mesh_initial_left
    twist_cp = [-3.75000, -0.943600, 6.71660]


# Create a dictionary with info and options about the aerodynamic
# lifting surface
surface = {
    # Wing definition
    "name": "wing",  # name of the surface
    "symmetry": True,  # if true, model one half of wing
    # reflected across the plane y = 0
    "S_ref_type": "wetted",  # how we compute the wing area,
    # can be 'wetted' or 'projected'
    "fem_model_type": "tube",
    "twist_cp": twist_cp,
    "mesh": mesh,
    # Aerodynamic performance of the lifting surface at
    # an angle of attack of 0 (alpha=0).
    # These CL0 and CD0 values are added to the CL and CD
    # obtained from aerodynamic analysis of the surface to get
    # the total CL and CD.
    # These CL0 and CD0 values do not vary wrt alpha.
    "CL0": 0.0,  # CL of the surface at alpha=0
    "CD0": 0.015,  # CD of the surface at alpha=0
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar
    # flow, used for viscous drag
    "t_over_c_cp": np.array([0.15]),  # thickness over chord ratio (NACA0015)
    "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
    # thickness
    "with_viscous": True,  # if true, compute viscous drag
    "with_wave": True,  # if true, compute wave drag
}

# Create the OpenMDAO problem
prob = om.Problem()

# Create an independent variable component that will supply the flow
# conditions to the problem.
indep_var_comp = om.IndepVarComp()
indep_var_comp.add_output("v", val=248.136, units="m/s")
indep_var_comp.add_output("alpha", val=5.0, units="deg")
indep_var_comp.add_output("Mach_number", val=0.84)
indep_var_comp.add_output("re", val=1.0e6, units="1/m")
indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")

# Add this IndepVarComp to the problem model
prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

# Create and add a group that handles the geometry for the
# aerodynamic lifting surface
geom_group = Geometry(surface=surface)
prob.model.add_subsystem(surface["name"], geom_group)

# Create the aero point group, which contains the actual aerodynamic
# analyses
aero_group = AeroPoint(surfaces=[surface])
point_name = "aero_point_0"
prob.model.add_subsystem(point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])

name = surface["name"]

# Connect the mesh from the geometry component to the analysis point
prob.model.connect(name + ".mesh", point_name + "." + name + ".def_mesh")

# Perform the connections with the modified names within the
# 'aero_states' group.
prob.model.connect(name + ".mesh", point_name + ".aero_states." + name + "_def_mesh")

prob.model.connect(name + ".t_over_c", point_name + "." + name + "_perf." + "t_over_c")

# Import the Scipy Optimizer and set the driver of the problem to use
# it, which defaults to an SLSQP optimization method
import openmdao.api as om

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["tol"] = 1e-9

recorder = om.SqliteRecorder("aero_analysis_test.db")
prob.driver.add_recorder(recorder)
prob.driver.recording_options["record_derivatives"] = True
prob.driver.recording_options["includes"] = ["*"]

# Setup problem and add design variables, constraint, and objective
prob.model.add_design_var("wing.twist_cp", lower=-10.0, upper=15.0)
prob.model.add_constraint(point_name + ".wing_perf.CL", equals=0.5)
prob.model.add_objective(point_name + ".wing_perf.CD", scaler=1e4)

# Set up and run the optimization problem
prob.setup()
prob.run_model()
# n2(prob)






###########################################################################
#######                    Comparison Plot                          #######
###########################################################################


WingMeshPlot = plt.figure(1)
ax = plt.axes(projection='3d')

# plot my mesh
nx = mesh_initial_left.shape[0]
ny_total = mesh_initial_left.shape[1]

# plot intial mesh (no twist)
# spanwise lines from the leading edge
for i in range(0, nx):
    x_spanwise = mesh_initial_left[i,:,0]
    y_spanwise = mesh_initial_left[i,:,1]
    z_spanwise = mesh_initial_left[i,:,2]
    ax.plot(x_spanwise, y_spanwise, z_spanwise, color='blue')
# chordwise lines from the tip 
for i in range(0, ny_total):
    x_chordwise = mesh_initial_left[:,i,0]
    y_chordwise = mesh_initial_left[:,i,1]
    z_chordwise = mesh_initial_left[:,i,2]
    ax.plot(x_chordwise, y_chordwise, z_chordwise, color='blue')

# plot def_mesh (actual meshed used in analysis)
# spanwise lines from the leading edge
for i in range(0, nx):
    x_spanwise = mesh_deformed_left[i,:,0]
    y_spanwise = mesh_deformed_left[i,:,1]
    z_spanwise = mesh_deformed_left[i,:,2]
    ax.plot(x_spanwise, y_spanwise, z_spanwise, color='cyan')
# chordwise lines from the tip 
for i in range(0, ny_total):
    x_chordwise = mesh_deformed_left[:,i,0]
    y_chordwise = mesh_deformed_left[:,i,1]
    z_chordwise = mesh_deformed_left[:,i,2]
    ax.plot(x_chordwise, y_chordwise, z_chordwise, color='cyan')



###########################################################################
# plot mesh generated by OpenAeroStruct
mesh_initial_left_original = prob["wing.mesh.taper.mesh"]
if option == 0:
    shift = mesh_initial_left_original[0,-1,0]
else:
    shift = 0 #in this case original OpenAeroStruct also uses my mesh, therefore no need to shift

mesh_initial_left_original[:,:,0] = mesh_initial_left_original[:,:,0] - shift # shift to same location
nx = mesh_initial_left_original.shape[0]
ny_total = mesh_initial_left_original.shape[1]

# plot intial mesh (no twist)
# spanwise lines from the leading edge
for i in range(0, nx):
    x_spanwise = mesh_initial_left_original[i,:,0]
    y_spanwise = mesh_initial_left_original[i,:,1]
    z_spanwise = mesh_initial_left_original[i,:,2]
    ax.plot(x_spanwise, y_spanwise, z_spanwise, color='red')
# chordwise lines from the tip 
for i in range(0, ny_total):
    x_chordwise = mesh_initial_left_original[:,i,0]
    y_chordwise = mesh_initial_left_original[:,i,1]
    z_chordwise = mesh_initial_left_original[:,i,2]
    ax.plot(x_chordwise, y_chordwise, z_chordwise, color='red')

mesh_deformed_left_original = prob["aero_point_0.wing.def_mesh"]
mesh_deformed_left_original[:,:,0] = mesh_deformed_left_original[:,:,0] - shift # shift to same location
nx = mesh_deformed_left_original.shape[0]
ny_total = mesh_deformed_left_original.shape[1]
# plot def_mesh (actual meshed used in analysis)
# spanwise lines from the leading edge
for i in range(0, nx):
    x_spanwise = mesh_deformed_left_original[i,:,0]
    y_spanwise = mesh_deformed_left_original[i,:,1]
    z_spanwise = mesh_deformed_left_original[i,:,2]
    ax.plot(x_spanwise, y_spanwise, z_spanwise, color='pink')
# chordwise lines from the tip 
for i in range(0, ny_total):
    x_chordwise = mesh_deformed_left_original[:,i,0]
    y_chordwise = mesh_deformed_left_original[:,i,1]
    z_chordwise = mesh_deformed_left_original[:,i,2]
    ax.plot(x_chordwise, y_chordwise, z_chordwise, color='pink')

plt.show()



###########################################################################
#######                    Comparison Print                         #######
###########################################################################

fileName = "Results/Model_Execution/SingleRunAero_Validation.txt"
f = open(fileName, "w")

mesh_initial_left_str = np.array2string(mesh_initial_left, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
f.write("\nmy mesh_initial_left   = ")
f.write(mesh_initial_left_str + "\n")

mesh_initial_left_original_str = np.array2string(mesh_initial_left_original, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
f.write("\nOSA mesh_initial_left  = ")
f.write(mesh_initial_left_original_str + "\n")


mesh_deformed_left_str = np.array2string(mesh_deformed_left, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
f.write("\nmy mesh_deformed_left  = ")
f.write(mesh_deformed_left_str + "\n")

mesh_deformed_left_original_str = np.array2string(mesh_deformed_left_original, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
f.write("\nOSA mesh_deformed_left = ")
f.write(mesh_deformed_left_original_str + "\n")


sec_forces_str = np.array2string(sec_forces, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
f.write("\nmy sec_forces  = ")
f.write(sec_forces_str + "\n")

sec_forces_original = prob["aero_point_0.aero_states.wing_sec_forces"]
sec_forces_original_str = np.array2string(sec_forces_original, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
f.write("\nOSA sec_forces = ")
f.write(sec_forces_original_str + "\n")


f.write("\nmy CL          = ")
f.write("{:.6f}".format(CL_Wing_total))

CL_original = prob["aero_point_0.CL"][0]
f.write("\nOSA CL         = ")
f.write("{:.6f}".format(CL_original))


f.write("\nmy CD          = ")
f.write("{:.6f}".format(CD_Wing_total))

CD_original = prob["aero_point_0.CD"][0]
f.write("\nOSA CD         = ")
f.write("{:.6f}".format(CD_original))


f.write("\nmy CDw         = ")
f.write("{:.6f}".format(CD_Wing_w))

CDw_original = prob["aero_point_0.wing_perf.CDw"][0]
f.write("\nOSA CDw        = ")
f.write("{:.6f}".format(CDw_original))


f.write("\nmy L           = ")
f.write("{:.6f}".format(L_Wing))

L_Wing_Original = prob["aero_point_0.wing_perf.L"][0]
f.write("\nOSA L          = ")
f.write("{:.6f}".format(L_Wing_Original))


f.write("\nmy D           = ")
f.write("{:.6f}".format(D_Wing))

D_Wing_Original = prob["aero_point_0.wing_perf.D"][0]
f.write("\nOSA D          = ")
f.write("{:.6f}".format(D_Wing_Original))

print("done")
