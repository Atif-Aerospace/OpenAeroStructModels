"""
Runs the actual aerodynamics analysis using OpenAeroStruct AeroPoint component
"""

from openaerostruct.aerodynamics.aero_groups import AeroPoint
from openaerostruct.transfer.load_transfer import LoadTransfer
from openmdao.api import n2
import openmdao.api as om
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def AerodynamicsAnalysis(wingInfo, velocity, alpha, Mach_number, Re, rho, cg_location, mesh_deformed_left, t_over_c):

    # Create the OpenMDAO problem
    prob = om.Problem()

    # Create an independent variable component that will supply the flow
    # conditions to the problem.
    flight_conditions_comp = om.IndepVarComp()
    flight_conditions_comp.add_output("v", val = velocity, units="m/s")
    flight_conditions_comp.add_output("alpha", val = alpha, units="deg")
    flight_conditions_comp.add_output("Mach_number", val = Mach_number)
    flight_conditions_comp.add_output("re", val = Re, units="1/m")
    flight_conditions_comp.add_output("rho", val = rho, units="kg/m**3")
    flight_conditions_comp.add_output("cg", val = cg_location, units="m")

    # Add this IndepVarComp to the problem model
    prob.model.add_subsystem("flight_conditions", flight_conditions_comp, promotes=["*"])

    # No need for this Geometry group which created in WingDefinition.py
    # winggeom_group = Geometry(surface=surf_dict)
    # prob.model.add_subsystem("winggeom_group", winggeom_group)

    # Create an independent variable component to receive 
    # def_mesh and t_over_c which are produced in WingDefinition.py
    input_var_comp = om.IndepVarComp()
    # the values of def_mesh should be the same as from structure design
    input_var_comp.add_output("def_mesh", val = mesh_deformed_left, units="m") 

    input_var_comp.add_output("t_over_c", val = t_over_c)

    # Add this IndepVarComp to the problem model
    prob.model.add_subsystem("input_var", input_var_comp)

    # Create the aero point group, which contains the actual aerodynamic
    # analyses
    aeropoint_group = AeroPoint(surfaces=[wingInfo])
    prob.model.add_subsystem("aeropoint_group", aeropoint_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])

    # Connect the mesh from the geometry component to the analysis point
    prob.model.connect("input_var.def_mesh", "aeropoint_group.wing.def_mesh")

    # Perform the connections with the modified names within the
    # 'aero_states' group.
    prob.model.connect("input_var.def_mesh", "aeropoint_group.aero_states.wing_def_mesh")
    prob.model.connect("input_var.t_over_c", "aeropoint_group.wing_perf.t_over_c")

    # Add a loads component to the coupled group
    prob.model.add_subsystem("loadtransfer_group", LoadTransfer(surface=wingInfo))
    prob.model.connect("input_var.def_mesh", "loadtransfer_group.def_mesh")
    prob.model.connect("aeropoint_group" + ".aero_states.wing_sec_forces", "loadtransfer_group.sec_forces")

    # Set up and run the model
    prob.setup()
    prob.run_model()
    # n2(prob)

    # =================================== Variables of interest ===================================
    # Wing reference area
    wing_area_asref = prob["aeropoint_group.wing.S_ref"][0]

    # Lift and Drag Forces of the Wing
    L_Wing = prob["aeropoint_group.wing_perf.L"][0]
    D_Wing = prob["aeropoint_group.wing_perf.D"][0]

    # Lift and Drag Coefficients of the wing
    CL_Wing_total = prob["aeropoint_group.wing_perf.CL"][0] # induced drag
    CD_Wing_i = prob["aeropoint_group.wing_perf.CDi"][0] # induced drag
    CD_Wing_v = prob["aeropoint_group.wing_perf.CDv"][0] # viscous drag
    CD_Wing_w = prob["aeropoint_group.wing_perf.CDw"][0] # wave drag
    CD_Wing_total = prob["aeropoint_group.wing_perf.CD"][0] # total drag

    # Moments Coefficients
    CM_Wing = np.copy(prob["aeropoint_group.CM"])
    CM_Wing_roll = CM_Wing[0]
    CM_Wing_pitch = CM_Wing[1]
    CM_Wing_yaw = CM_Wing[2]

    # Lift over Drag
    LoD_Wing = CL_Wing_total/CD_Wing_total

    # Array containing the sectional forces acting on each panel.
    sec_forces = np.copy(prob["loadtransfer_group.sec_forces"])

    # break sec_forces into forces along and moments around x,y,z axies
    sec_forces_Fx = sec_forces[:,:,0]
    sec_forces_Fy = sec_forces[:,:,1]
    sec_forces_Fz = sec_forces[:,:,2]


    # loads applied on the FEM component at each node
    # The first 3 indices are N and the last 3 are N*m
    loads = np.copy(prob["loadtransfer_group.loads"])

    # From load_transfer.py
    # line 145: The moment arm is between the aerodynamic centers of each panel and the FEM nodes.
    # line 132: Compute the aerodynamic centers at the quarter-chord point of each panel 
    # line 141: The structural nodes is computed based on the fem_origin location (weighted sum of the LE and TE mesh vertices)

    # break loads into forces along and moments around x,y,z axies
    # these values are negative which are strange but in the provided example Fx are also negativte, see also 
    # From D:\Tools\OpenAeroStruct\openaerostruct\aerodynamics\lift_drag.py line 59 and 62
    # we can see that negative direction of x axis is pointing forward?
    loads_Fx = loads[:,0]
    loads_Fy = loads[:,1]
    loads_Fz = loads[:,2]
    loads_Mx = loads[:,3]
    loads_My = loads[:,4]
    loads_Mz = loads[:,5]

    # Temp
    # double the root load
    loads_Fx[-1] = loads_Fx[-1] * 2
    loads_Fy[-1] = loads_Fy[-1] * 2
    loads_Fz[-1] = loads_Fz[-1] * 2
    # Temp

    # Produce Drag Polar
    # Create arrays to hold polar data
    N_Mach = 10
    N_alpha = 10
    N_CL = 10

    Mach_DP = np.linspace(0.0, 0.9, N_Mach)
    alpha_DP = np.linspace(-5, 10, N_alpha)

    CL_DP_Raw = np.zeros((N_Mach, N_alpha))
    CD_total_DP_Raw = np.zeros((N_Mach, N_alpha))
    CD_i_DP_raw = np.zeros((N_Mach, N_alpha))

    # CL_DP_fixed = np.linspace(0.0, 0.50, N_CL)
    CL_DP_fixed = np.array([0.00000, 0.05000, 0.10000, 0.15000, 0.20000, 0.25000, 0.30000, 0.40000, 0.50000, 0.60000]) 
    CD_DP_final = np.zeros((N_Mach, N_CL))
    CD_i_DP_final = np.zeros((N_Mach, N_CL))

    # Sweep through alphas and Mach to create polar
    for i in range(len(Mach_DP)):

        for j in range(len(alpha_DP)):
            # Set alpha
            prob["Mach_number"] = Mach_DP[i] #check if v needs to be changed as well
            # prob["v"] = Mach_DP[i] * 296.54 #check if v needs to be changed as well
            prob["alpha"] = alpha_DP[j]


            # Run analysis
            prob.run_model()

            # Record CL, CD
            CL_DP_Raw[i,j] = prob["aeropoint_group.wing_perf.CL"][0]
            CD_total_DP_Raw[i,j] = prob["aeropoint_group.wing_perf.CD"][0]
            CD_i_DP_raw[i,j] = prob["aeropoint_group.wing_perf.CDi"][0] # induced drag

        # Using Interpolate to get CD at fixed CLs
        tck_total = interpolate.splrep(CL_DP_Raw[i,:], CD_total_DP_Raw[i,:], s=0)
        CD_DP_final[i,:] = interpolate.splev(CL_DP_fixed, tck_total, der=0)

        tck_induced = interpolate.splrep(CL_DP_Raw[i,:], CD_i_DP_raw[i,:], s=0)
        CD_i_DP_final[i,:] = interpolate.splev(CL_DP_fixed, tck_induced, der=0)

        # Plot polar
        plotSwitch = False
        if plotSwitch == True:
            plt.figure()
            plt.plot(CD_total_DP_Raw[i,:] * 1e4, CL_DP_Raw[i,:], 'x')
            plt.plot(CD_DP_final[i,:] * 1e4, CL_DP_fixed, "-o")
            plt.grid(color="lightgray", linestyle="-", linewidth=1)
            plt.legend(['True', 'Interpolate', 'True'])
            plt.xlabel("$C_D$ (counts)")
            plt.ylabel("$C_L$")
            plt.title("Drag polar")
            plt.show()
            # plt.show(block=False)

    return (wing_area_asref, L_Wing, D_Wing, LoD_Wing,
    CL_Wing_total, CD_Wing_i, CD_Wing_v, CD_Wing_w, CD_Wing_total, CD_Wing_total,
    CL_DP_fixed, CD_DP_final, CD_i_DP_final,
    CM_Wing, CM_Wing_roll, CM_Wing_pitch, CM_Wing_yaw,
    sec_forces, sec_forces_Fx, sec_forces_Fy, sec_forces_Fz,
    loads, loads_Fx, loads_Fy, loads_Fz, loads_Mx, loads_My, loads_Mz)