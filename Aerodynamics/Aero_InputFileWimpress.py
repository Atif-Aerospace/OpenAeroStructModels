################### Planform Definition ###################

wing_area_wpref             = 383.740000    # [m^2]     Reference area based on Wimpress defintion
aspect_ratio_wpref          =   9.000000    # [ ]       wing aspect ratio
kink_location_ratio         =   0.370000    # [ ]       spanwise location of the kink from the center line / semispan
body_side_ratio             =   0.100000    # [ ]       spanwise location of the side of fuselage from the center line / semispan
taper_ratio_trap            =   0.275000    # [ ]       Taper ratio at tip/root for the trap-wing area
root_chord_extension_ratio  =   1.372874    # [ ]       Actrual root chord / trap root chord
trap_quarter_sweep          =  35.000000    # [deg]     Trapezoidal area leading-edge sweep angle


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

################### Assumed variables ###################

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

###################  Assumed Deformation ###################

# Assumed deformation on leading edge
# index 1, 2, 3 are from tip to root
delta_x_LE_1        = 0.009000
delta_x_LE_2        = 0.004000
delta_x_LE_3        = 0.001000  

delta_y_LE_1        = 0.000090
delta_y_LE_2        = 0.000040
delta_y_LE_3        = 0.000010

delta_z_LE_1        = 0.090000
delta_z_LE_2        = 0.040000
delta_z_LE_3        = 0.010000

# Assumed deformation on trailing edge
# index 1, 2, 3 are from tip to root
delta_x_TE_1        = 0.009000
delta_x_TE_2        = 0.004000
delta_x_TE_3        = 0.001000  

delta_y_TE_1        = 0.000090
delta_y_TE_2        = 0.000040
delta_y_TE_3        = 0.000010

delta_z_TE_1        = 0.090000
delta_z_TE_2        = 0.040000
delta_z_TE_3        = 0.010000



