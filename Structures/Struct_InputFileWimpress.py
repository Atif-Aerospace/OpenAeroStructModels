import numpy as np
################### Material Properties ###################

E                   = 70.0e9            # [Pa] Young's modulus of the spar
G                   = 30.0e9            # [Pa] shear modulus of the spar
yieldStress         = 500.0e6 / 2.5     # [Pa] yield stress divided by 2.5 for limiting case
mrho                = 3.0e3             # [kg/m^3] material density



################### Assumed variables ###################

wing_weight_ratio                 =   2.000000    # [ ]      

################### Planform Definition ###################

wing_area_wpref             = 383.740000    # [m^2]     wing area in the trapezoid area in the Wimpress definition
aspect_ratio_wpref          =   9.000000    # [ ]       wing aspect ratio
kink_location_ratio         =   0.370000    # [ ]       spanwise location of the kink from the center line / semispan
body_side_ratio             =   0.100000    # [ ]       spanwise location of the side of fuselage from the center line / semispan
taper_ratio_trap            =   0.275000    # [ ]       Taper ratio at tip/root for the trap-wing area
root_chord_extension_ratio  =   1.372874    # [ ]       Actrual root chord / trap root chord
trap_quarter_sweep          =  35.000000    # [deg]     inboard leading-edge sweep angle


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

# chordwise location of maximum (NACA0015) thickness
c_max_t             =   0.303000    # [ ]       

############ Structure Design Variables ############

fem_origin          =   0.350000    # [ ]       normalized chordwise location of the spar
thickness_cp_1      =   0.05
thickness_cp_2      =   0.05
thickness_cp_3      =   0.05

####################### Loads ######################
# Loads on the structural (spar) nodes, from tip to root:
#  Forces in [N]
#  Bending Moments in [N.m]
# loads_Fx            =  [-3048.874830, -8049.131076, -11794.825236, -14964.573007, -17386.897108, -14545.454519, -10657.122856]
# loads_Fy            =  [3226.235167, 7365.270459, 8407.813943, 7916.857034, 7311.487887, 6269.430025, 5212.041039]
# loads_Fz            =  [28543.906354, 76509.644253, 116952.174321, 160342.239542, 229686.530966, 300342.680528, 324023.905362]
# loads_Mx            =  [67474.346278, 45053.664144, 48362.409955, 51369.335294, 164944.105160, 65472.998931, -439781.937673]
# loads_My            =  [56521.526875, 61748.554484, 92414.612971, 133329.278205, 242572.065277, 319553.435956, -76193.955789]
# loads_Mz            =  [1620.836443, 3241.671218, 6658.228366, 12134.503883, 24877.420696, 24477.362453, 10003.414330]

loads_Fx = [-1029.064316, -2531.531482, -3403.647142, -4180.870113, -4921.170350, -5620.316368, -6257.801360, -6802.043976, -7206.724223, -7541.507961, -7488.627829, -6755.896025, -5344.778160, -2813.191846, -1187.088780]
loads_Fy = [1109.663594, 2566.086938, 3109.517540, 3423.229171, 3588.135955, 3616.061850, 3510.741859, 3282.945657, 2962.360903, 2783.819639, 2613.287960, 2291.244142, 2122.827000, 2514.957894, 2923.620965]
loads_Fz = [9384.136225, 23132.306068, 31326.918304, 39009.790113, 46827.612275, 54871.156340, 63115.333452, 71512.036577, 80056.200514, 91967.271453, 104735.916590, 114288.414243, 121935.488654, 127065.513386, 128929.764333]
loads_Mx = [9881.086668, 4540.904096, 3943.123512, 3934.728773, 4028.089519, 4130.754578, 4224.138767, 4316.129688, 4458.275275, 10886.208770, 5899.050473, 5251.319562, 4285.271093, 2970.076305, -69751.722740]
loads_My = [10298.314341, 11649.167458, 15069.256222, 19856.140781, 25682.446100, 32468.023206, 40145.843583, 48562.679764, 57305.881758, 68164.207532, 89862.219594, 109952.734979, 124518.469145, 123619.878735, 19459.287606]
loads_Mz = [100.239061, 19.178273, 318.701527, 753.774651, 1379.543610, 2236.344417, 3351.026857, 4733.413377, 6378.021573, 8883.844648, 11027.077047, 13661.473375, 15678.618751, 15288.175853, 8682.661131]


ny = len(loads_Fx)
loads = np.ones((ny, 6)) * 1e5
loads[:,0] = loads_Fx
loads[:,1] = loads_Fy
loads[:,2] = loads_Fz
loads[:,3] = loads_Mx
loads[:,4] = loads_My
loads[:,5] = loads_Mz


