from modelsCloudDeployment_V3 import OAS_Aerodynamics_V3 as Aero_V3
from modelsCloudDeployment_V3 import OAS_Structure_V3 as Struct_V3
from modelsCloudDeployment_V3 import OAS_Iteration_V3 as IAS_V3
import numpy as np



################### Planform Definition ###################

WingArea                    = 383.740000    # [m^2]     Reference area based on Wimpress defintion
WingAspectRatio             =   9.000000    # [ ]       wing aspect ratio
kink_location_ratio         =   0.370000    # [ ]       spanwise location of the kink from the center line / semispan
body_side_ratio             =   0.100000    # [ ]       spanwise location of the side of fuselage from the center line / semispan
WingTaperRatio              =   0.275000    # [ ]       Taper ratio at tip/root for the trap-wing area
root_chord_extension_ratio  =   1.372874    # [ ]       Actrual root chord / trap root chord
WingSweep25                 =  35.000000    # [deg]     Trapezoidal area leading-edge sweep angle


################### 3D geometry ###################

dihedral            =   1.500000    # [deg]     wing dihedral

# Control points of twist and thickness to chord ratio. 
# index 1, 2, 3 are from TIP to ROOT
twist_cp_1          =  -3.000000    # [deg]     
twist_cp_2          =   0.000000    # [deg]     
twist_cp_3          =   3.000000    # [deg]
t_over_c_cp_1       =   0.150000    # [ ]
t_over_c_cp_2       =   0.150000    # [ ]
t_over_c_cp_3       =   0.150000    # [ ]

################### Aero Assumed variables ###################

CL0                 =   0.000000    # [ ]       CL of the surface at alpha=0
CD0                 =   0.015000    # [ ]       CD of the surface at alpha=0
k_lam               =   0.050000    # [ ]       percentage of chord with laminar flow, used for viscous drag 
c_max_t             =   0.303000    # [ ]       chordwise location of maximum (NACA0015) thickness
fem_origin          =   0.350000    # [ ]       normalized chordwise location of the spar
cg_location_x       =   0.000000    # [m]       centre of gravity location in x axis
cg_location_y       =   0.000000    # [m]       centre of gravity location in y axis
cg_location_z       =   0.000000    # [m]       centre of gravity location in z axis

################### Flight Conditions ###################

velocity            = 248.136000    # [m/s]     flight speed   sonic speed 296.54 m /s
alpha               =   5.000000    # [deg]     Angle of attack
Mach_number         =   0.840000    # [ ]       Mach Number
Re                  =   1000000     # [1/m]     Unit? Reynolds Number 
rho                 =   0.380000    # [kg/m**3] Air density   35000ft


################### Material Properties ###################

E                   = 70.0e9            # [Pa] Young's modulus of the spar
G                   = 30.0e9            # [Pa] shear modulus of the spar
yieldStress         = 500.0e6 / 2.5     # [Pa] yield stress divided by 2.5 for limiting case
mrho                = 3.0e3             # [kg/m^3] material density


############ Structure Design Variables ############
thickness_cp_1      =   0.05            # [m] Control points of spar thickness TIP
thickness_cp_2      =   0.05            # [m] Control points of spar thickness MIDDLE
thickness_cp_3      =   0.05            # [m] Control points of spar thickness ROOT

################### Structure Assumed variables ###################

wing_weight_ratio                 =   2.000000    # [ ]      

dragPolar = 1

(LoD_Wing, CL_Wing_total, dragPolarCL, dragPolarCD, CD_i_dragPolar_intp, 
    thickness_intersects, element_mass, wingWeight, vonmisesStress, failure, # If we need to plot uncomment the lines below
    loads_Converged, 
    mesh_undeformed_left_x_Converged, mesh_undeformed_left_y_Converged, mesh_undeformed_left_z_Converged,
    mesh_deformed_left_x_Converged, mesh_deformed_left_y_Converged, mesh_deformed_left_z_Converged,
    mesh_delta_left_x_Converged, mesh_delta_left_y_Converged, mesh_delta_left_z_Converged, 
    spar_pts_undeformed_left_x_Converged, spar_pts_undeformed_left_y_Converged, spar_pts_undeformed_left_z_Converged,
    spar_pts_deformed_left_x_Converged, spar_pts_deformed_left_y_Converged, spar_pts_deformed_left_z_Converged) = IAS_V3(WingArea, WingAspectRatio, kink_location_ratio, body_side_ratio,   # Planform definition  
        WingTaperRatio, root_chord_extension_ratio, WingSweep25,   # Planform definition
        dihedral, twist_cp_1, twist_cp_2, twist_cp_3, t_over_c_cp_1, t_over_c_cp_2, t_over_c_cp_3,  # 3D geometry: 1, 2, 3 are from tip to root
        CL0, CD0, k_lam, c_max_t, fem_origin,   # Assumed variables
        velocity, alpha, Mach_number, Re, rho,  # Flight conditions
        cg_location_x, cg_location_y, cg_location_z,    # cg locations (if known)
        thickness_cp_1, thickness_cp_2, thickness_cp_3,
        E, G, yieldStress, mrho, wing_weight_ratio, dragPolar)


print("done")



