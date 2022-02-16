"""train seperate ANN for different outputs"""
import numpy as np
from smt.sampling_methods import LHS

import pandas as pd
import matplotlib.pyplot as plt
import csv

import Aerodynamics.Aero_InputFileWimpress as aif # Read Input File
import Aerodynamics.Aero_IntegratedWorkflow as aiw
import Aerodynamics.Aero_InputFileWimpress as aif # Read Input File
import Auxiliary.Aux_PlotWingAndLoads as xpw
import Auxiliary.Aux_AeroWriteToFile as xwf
import NeuralNetworks.TrainNeuralNetworks as tnn



# constants
wing_area_wpref     = 383.740000    # [m^2]     wing area in the trapezoid area in the Wimpress definition
body_side_ratio     =   0.100000    # [ ]       spanwise location of the side of fuselage from the center line / semispan
CL0                 =   0.000000    # [ ]       CL of the surface at alpha=0
CD0                 =   0.015000    # [ ]       CD of the surface at alpha=0
k_lam               =   0.050000    # [ ]       percentage of chord with laminar flow, used for viscous drag 
c_max_t             =   0.303000    # [ ]       chordwise location of maximum (NACA0015) thickness
fem_origin          =   0.350000    # [ ]       normalized chordwise location of the spar
cg_location_x       =   0.000000    # [m]       centre of gravity location in x axis
cg_location_y       =   0.000000    # [m]       centre of gravity location in y axis
cg_location_z       =   0.000000    # [m]       centre of gravity location in z axis
velocity            = 248.136000    # [m/s]     flight speed
alpha               =   5.000000    # [deg]     Angle of attack
Mach_number         =   0.840000    # [ ]       Mach Number
Re                  =   1000000     # [1/m]     Unit? Reynolds Number 
rho                 =   0.380000    # [kg/m**3] Air density

# initiate the upper and lower bound for design of experiment
xlimits = np.array([
[7.0, 10.0],     # [ ]   aspect_ratio_wpref:         wing aspect ratio
[0.37, 0.40],    # [ ]   kink_location_ratio:        spanwise location of the kink from the center line / semispan
[0.25, 0.30],    # [ ]   taper_ratio_trap:           Taper ratio at tip/root for the trap-wing area
[1.20, 1.50],    # [ ]   root_chord_extension_ratio: Actrual root chord / trap root chord
[33.00, 37.0],   # [deg] trap_quarter_sweep:         inboard leading-edge sweep angle

[0.00, 3.00],    # [deg] dihedral:                   wing dihedral
[-5.00, 5.00],   # [deg] twist_cp_1:                 Control points of twist and thickness to chord ratio. index 1, 2, 3 are from TIP to ROOT 
[-5.00, 5.00],   # [deg] twist_cp_2                       
[-5.00, 5.00],   # [deg] twist_cp_3 

[0.10, 0.15],    # [ ]   t_over_c_cp_1
[0.10, 0.15],    # [ ]   t_over_c_cp_2
[0.10, 0.15],    # [ ]   t_over_c_cp_3

[-0.05, 0.05],   # [ ]   delta_x_LE_1:               Assumed deformation on leading edge. index 1, 2, 3 are from tip to root
[-0.05, 0.05],   # [ ]   delta_x_LE_2
[-0.05, 0.05],   # [ ]   delta_x_LE_3

[-0.05, 0.05],   # [ ]   delta_y_LE_1
[-0.05, 0.05],   # [ ]   delta_y_LE_2
[-0.05, 0.05],   # [ ]   delta_y_LE_3

[-0.05, 0.05],   # [ ]   delta_z_LE_1
[-0.05, 0.05],   # [ ]   delta_z_LE_2
[-0.05, 0.05],   # [ ]   delta_z_LE_3

[-0.05, 0.05],   # [ ]   delta_x_TE_1:               Assumed deformation on trailing edge, index 1, 2, 3 are from tip to root
[-0.05, 0.05],   # [ ]   delta_x_TE_2
[-0.05, 0.05],   # [ ]   delta_x_TE_3

[-0.05, 0.05],   # [ ]   delta_y_TE_1
[-0.05, 0.05],   # [ ]   delta_y_TE_2
[-0.05, 0.05],   # [ ]   delta_y_TE_3

[-0.05, 0.05],   # [ ]   delta_z_TE_1
[-0.05, 0.05],   # [ ]   delta_z_TE_2
[-0.05, 0.05],  # [ ]   delta_z_TE_3
    ])

# Produce normalized samples
numTrain = 6400
sampling = LHS(xlimits=xlimits)
inputTrain = sampling(numTrain)
print(inputTrain.shape)

# Planform definition 
aspect_ratio_wpref          = inputTrain[:,0]
kink_location_ratio         = inputTrain[:,1]
taper_ratio_trap            = inputTrain[:,2]
root_chord_extension_ratio  = inputTrain[:,3]
trap_quarter_sweep          = inputTrain[:,4]

# 3D geometry: 1, 2, 3 are from tip to root
dihedral                    = inputTrain[:,5]
twist_cp_1                  = inputTrain[:,6]
twist_cp_2                  = inputTrain[:,7]
twist_cp_3                  = inputTrain[:,8]
t_over_c_cp_1               = inputTrain[:,9]
t_over_c_cp_2               = inputTrain[:,10]
t_over_c_cp_3               = inputTrain[:,11]

# Assumed deformation on leading edge: 1, 2, 3 are from tip to root
delta_x_LE_1                = inputTrain[:,12]
delta_x_LE_2                = inputTrain[:,13]
delta_x_LE_3                = inputTrain[:,14]
delta_y_LE_1                = inputTrain[:,15]
delta_y_LE_2                = inputTrain[:,16]
delta_y_LE_3                = inputTrain[:,17]
delta_z_LE_1                = inputTrain[:,18]
delta_z_LE_2                = inputTrain[:,19]
delta_z_LE_3                = inputTrain[:,20]

# Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
delta_x_TE_1                = inputTrain[:,21]
delta_x_TE_2                = inputTrain[:,22]
delta_x_TE_3                = inputTrain[:,23]
delta_y_TE_1                = inputTrain[:,24]
delta_y_TE_2                = inputTrain[:,25]
delta_y_TE_3                = inputTrain[:,26]
delta_z_TE_1                = inputTrain[:,27]
delta_z_TE_2                = inputTrain[:,28]
delta_z_TE_3                = inputTrain[:,29]


# Record samples
outputTrain = np.zeros([numTrain, 44])


#  Run Actual aerodynamic analysis
for i in range(numTrain):
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
        ) = aiw.AerodynamicsFunction_wp_10_CP(wing_area_wpref, aspect_ratio_wpref[i], kink_location_ratio[i], body_side_ratio,   # Planform definition  
        taper_ratio_trap[i], root_chord_extension_ratio[i], trap_quarter_sweep[i],   # Planform definition
        dihedral[i], twist_cp_1[i], twist_cp_2[i], twist_cp_3[i], t_over_c_cp_1[i], t_over_c_cp_2[i], t_over_c_cp_3[i],  # 3D geometry: 1, 2, 3 are from tip to root
        CL0, CD0, k_lam, c_max_t, fem_origin,   # Assumed variables
        velocity, alpha, Mach_number, Re, rho,  # Flight conditions
        cg_location_x, cg_location_y, cg_location_z,    # cg locations (if known)
        delta_x_LE_1[i], delta_x_LE_2[i], delta_x_LE_3[i],   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
        delta_y_LE_1[i], delta_y_LE_2[i], delta_y_LE_3[i],
        delta_z_LE_1[i], delta_z_LE_2[i], delta_z_LE_3[i],
        delta_x_TE_1[i], delta_x_TE_2[i], delta_x_TE_3[i],   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
        delta_y_TE_1[i], delta_y_TE_2[i], delta_y_TE_3[i],
        delta_z_TE_1[i], delta_z_TE_2[i], delta_z_TE_3[i])

    outputTrain[i,0] = LoD_Wing
    outputTrain[i,1] = CL_Wing_total
    # mesh_deformed_left_x
    # mesh_deformed_left_y
    # mesh_deformed_left_z
    outputTrain[i,2:(2+ny_total)]                 = loads_Fx.flatten()
    outputTrain[i,(2+ny_total):(2+2*ny_total)]    = loads_Fy.flatten()
    outputTrain[i,(2+2*ny_total):(2+3*ny_total)]  = loads_Fz.flatten()
    outputTrain[i,(2+3*ny_total):(2+4*ny_total)]  = loads_Mx.flatten()
    outputTrain[i,(2+4*ny_total):(2+5*ny_total)]  = loads_My.flatten()
    outputTrain[i,(2+5*ny_total):(2+6*ny_total)]  = loads_Mz.flatten()
    
# Write training to files
rawDataTrain = np.concatenate((inputTrain, outputTrain), axis=1)
header = ['aspect_ratio_wpref', 'kink_location_ratio',
        'taper_ratio_trap', 'root_chord_extension_ratio', 'trap_quarter_sweep',
        'dihedral', 'twist_cp_1', 'twist_cp_2', 'twist_cp_3', 't_over_c_cp_1', 't_over_c_cp_2', 't_over_c_cp_3',
        'delta_x_LE_1', 'delta_x_LE_2', 'delta_x_LE_3',
        'delta_y_LE_1', 'delta_y_LE_2', 'delta_y_LE_3',
        'delta_z_LE_1', 'delta_z_LE_2', 'delta_z_LE_3',
        'delta_x_TE_1', 'delta_x_TE_2', 'delta_x_TE_3',
        'delta_y_TE_1', 'delta_y_TE_2', 'delta_y_TE_3',
        'delta_z_TE_1', 'delta_z_TE_2', 'delta_z_TE_3',
        'LoD_Wing', 'CL_Wing_total',
        'loads_Fx_1', 'loads_Fx_2', 'loads_Fx_3', 'loads_Fx_4', 'loads_Fx_5', 'loads_Fx_6', 'loads_Fx_7', 
        'loads_Fy_1', 'loads_Fy_2', 'loads_Fy_3', 'loads_Fy_4', 'loads_Fy_5', 'loads_Fy_6', 'loads_Fy_7', 
        'loads_Fz_1', 'loads_Fz_2', 'loads_Fz_3', 'loads_Fz_4', 'loads_Fz_5', 'loads_Fz_6', 'loads_Fz_7', 
        'loads_Mx_1', 'loads_Mx_2', 'loads_Mx_3', 'loads_Mx_4', 'loads_Mx_5', 'loads_Mx_6', 'loads_Mx_7', 
        'loads_My_1', 'loads_My_2', 'loads_My_3', 'loads_My_4', 'loads_My_5', 'loads_My_6', 'loads_My_7', 
        'loads_Mz_1', 'loads_Mz_2', 'loads_Mz_3', 'loads_Mz_4', 'loads_Mz_5', 'loads_Mz_6', 'loads_Mz_7']

with open('D:\Tools\OpenAeroStruct\openaerostruct\examples\mytest\FunctionalFormatTest_5\Results\ANN_Testing\TrainingSamples.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # write multiple rows
    writer.writerows(rawDataTrain)

batchSize = 320
epochNum1 = 3000
epochNum2 = 5000

Arch_LoD = np.array([30, 5])
myANN_LoD, OutMax_LoD, OutMin_LoD, OutDelta_LoD = tnn.TrainANN(Arch_LoD, inputTrain, outputTrain[:,0], 0, 'tanh', 0, batchSize, epochNum1)

Arch_CL = np.array([30, 5])
myANN_CL, OutMax_CL, OutMin_CL, OutDelta_CL = tnn.TrainANN(Arch_CL, inputTrain, outputTrain[:,1], 0, 'tanh', 0, batchSize, epochNum1)

Arch_Fx = np.array([45, 30, 15])
myANN_Fx, OutMax_Fx, OutMin_Fx, OutDelta_Fx = tnn.TrainANN(Arch_Fx, inputTrain, outputTrain[:,2:(2+ny_total)], 2, 'tanh', 0, batchSize, epochNum2)

Arch_Fy = np.array([30, 45, 30, 15])
myANN_Fy, OutMax_Fy, OutMin_Fy, OutDelta_Fy = tnn.TrainANN(Arch_Fy, inputTrain, outputTrain[:,(2+ny_total):(2+2*ny_total)], 2, 'tanh', 0, batchSize, epochNum2)

Arch_Fz = np.array([30, 45, 30, 15])
myANN_Fz, OutMax_Fz, OutMin_Fz, OutDelta_Fz = tnn.TrainANN(Arch_Fz, inputTrain, outputTrain[:,(2+2*ny_total):(2+3*ny_total)], 2, 'tanh', 0, batchSize, epochNum2)

Arch_Mx = np.array([30, 45, 30, 15])
myANN_Mx, OutMax_Mx, OutMin_Mx, OutDelta_Mx = tnn.TrainANN(Arch_Mx, inputTrain, outputTrain[:,(2+3*ny_total):(2+4*ny_total)], 2, 'tanh', 0, batchSize, epochNum2)

Arch_My = np.array([30, 45, 30, 15])
myANN_My, OutMax_My, OutMin_My, OutDelta_My = tnn.TrainANN(Arch_My, inputTrain, outputTrain[:,(2+4*ny_total):(2+5*ny_total)], 2, 'tanh', 0, batchSize, epochNum2)

Arch_Mz = np.array([30, 45, 30, 15])
myANN_Mz, OutMax_Mz, OutMin_Mz, OutDelta_Mz = tnn.TrainANN(Arch_Mz, inputTrain, outputTrain[:,(2+5*ny_total):(2+6*ny_total)], 2, 'tanh', 0, batchSize, epochNum2)



# Plot Test points 
numTest = 300
inputTest = sampling(numTest)
print(inputTest.shape)

# Planform definition 
aspect_ratio_wpref          = inputTest[:,0]
kink_location_ratio         = inputTest[:,1]
taper_ratio_trap            = inputTest[:,2]
root_chord_extension_ratio  = inputTest[:,3]
trap_quarter_sweep          = inputTest[:,4]

# 3D geometry: 1, 2, 3 are from tip to root
dihedral                    = inputTest[:,5]
twist_cp_1                  = inputTest[:,6]
twist_cp_2                  = inputTest[:,7]
twist_cp_3                  = inputTest[:,8]
t_over_c_cp_1               = inputTest[:,9]
t_over_c_cp_2               = inputTest[:,10]
t_over_c_cp_3               = inputTest[:,11]

# Assumed deformation on leading edge: 1, 2, 3 are from tip to root
delta_x_LE_1                = inputTest[:,12]
delta_x_LE_2                = inputTest[:,13]
delta_x_LE_3                = inputTest[:,14]
delta_y_LE_1                = inputTest[:,15]
delta_y_LE_2                = inputTest[:,16]
delta_y_LE_3                = inputTest[:,17]
delta_z_LE_1                = inputTest[:,18]
delta_z_LE_2                = inputTest[:,19]
delta_z_LE_3                = inputTest[:,20]

# Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
delta_x_TE_1                = inputTest[:,21]
delta_x_TE_2                = inputTest[:,22]
delta_x_TE_3                = inputTest[:,23]
delta_y_TE_1                = inputTest[:,24]
delta_y_TE_2                = inputTest[:,25]
delta_y_TE_3                = inputTest[:,26]
delta_z_TE_1                = inputTest[:,27]
delta_z_TE_2                = inputTest[:,28]
delta_z_TE_3                = inputTest[:,29]


# Record samples
outputTestActual = np.zeros([numTest, 44])


#  Run Actual aerodynamics analysis
for i in range(numTest):
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
        ) = aiw.AerodynamicsFunction_wp_10_CP(wing_area_wpref, aspect_ratio_wpref[i], kink_location_ratio[i], body_side_ratio,   # Planform definition  
        taper_ratio_trap[i], root_chord_extension_ratio[i], trap_quarter_sweep[i],   # Planform definition
        dihedral[i], twist_cp_1[i], twist_cp_2[i], twist_cp_3[i], t_over_c_cp_1[i], t_over_c_cp_2[i], t_over_c_cp_3[i],  # 3D geometry: 1, 2, 3 are from tip to root
        CL0, CD0, k_lam, c_max_t, fem_origin,   # Assumed variables
        velocity, alpha, Mach_number, Re, rho,  # Flight conditions
        cg_location_x, cg_location_y, cg_location_z,    # cg locations (if known)
        delta_x_LE_1[i], delta_x_LE_2[i], delta_x_LE_3[i],   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
        delta_y_LE_1[i], delta_y_LE_2[i], delta_y_LE_3[i],
        delta_z_LE_1[i], delta_z_LE_2[i], delta_z_LE_3[i],
        delta_x_TE_1[i], delta_x_TE_2[i], delta_x_TE_3[i],   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
        delta_y_TE_1[i], delta_y_TE_2[i], delta_y_TE_3[i],
        delta_z_TE_1[i], delta_z_TE_2[i], delta_z_TE_3[i])

    outputTestActual[i,0] = LoD_Wing
    outputTestActual[i,1] = CL_Wing_total
    # mesh_deformed_left_x
    # mesh_deformed_left_y
    # mesh_deformed_left_z
    outputTestActual[i,2:(2+ny_total)]                 = loads_Fx.flatten()
    outputTestActual[i,(2+ny_total):(2+2*ny_total)]    = loads_Fy.flatten()
    outputTestActual[i,(2+2*ny_total):(2+3*ny_total)]  = loads_Fz.flatten()
    outputTestActual[i,(2+3*ny_total):(2+4*ny_total)]  = loads_Mx.flatten()
    outputTestActual[i,(2+4*ny_total):(2+5*ny_total)]  = loads_My.flatten()
    outputTestActual[i,(2+5*ny_total):(2+6*ny_total)]  = loads_Mz.flatten()

# Write test samples to files
rawDataTest = np.concatenate((inputTest, outputTestActual), axis=1)
with open('D:\Tools\OpenAeroStruct\openaerostruct\examples\mytest\FunctionalFormatTest_5\Results\ANN_Testing\TestSamples.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # write multiple rows
    writer.writerows(rawDataTest)

# Predictions from each ANN
LoDPred = myANN_LoD.predict(inputTest) #- (-1)) / 2 * OutDelta_LoD + OutMin_LoD
CLPred  = myANN_CL.predict(inputTest ) #- (-1)) / 2 * OutDelta_CL + OutMin_CL
FxPred  = (myANN_Fx.predict(inputTest ) - (-1)) / 2 * OutDelta_Fx + OutMin_Fx
FyPred  = (myANN_Fy.predict(inputTest ) - (-1)) / 2 * OutDelta_Fy + OutMin_Fy
FzPred  = (myANN_Fz.predict(inputTest ) - (-1)) / 2 * OutDelta_Fz + OutMin_Fz
MxPred  = (myANN_Mx.predict(inputTest ) - (-1)) / 2 * OutDelta_Mx + OutMin_Mx
MyPred  = (myANN_My.predict(inputTest ) - (-1)) / 2 * OutDelta_My + OutMin_My
MzPred  = (myANN_Mz.predict(inputTest ) - (-1)) / 2 * OutDelta_Mz + OutMin_Mz
outputTestPred = np.concatenate((LoDPred, CLPred, FxPred, FyPred, FzPred, MxPred, MyPred, MzPred), axis=1)

# Write the predictions to file
rawDataTestPred = np.concatenate((inputTest, outputTestPred), axis=1)
with open('D:\Tools\OpenAeroStruct\openaerostruct\examples\mytest\FunctionalFormatTest_5\Results\ANN_Testing\TestSamplesPred.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # write multiple rows
    writer.writerows(rawDataTestPred)



# Write the error to file
outputTestError = (outputTestPred - outputTestActual) / outputTestActual
totalError = np.sum(abs(outputTestError),axis=1)
totalError = np.atleast_2d(totalError).T
rawDataTestCombine = np.concatenate((inputTest, outputTestActual, outputTestPred, outputTestError, totalError), axis=1)

header_pred = header[30:]
header_error = header[30:]
header_total = ['Total_Error']
for i in range(len(header_pred)):
    header_pred[i] = header_pred[i] + '_pred'
    header_error[i] = header_error[i] + '_error'

header_combined = header + header_pred + header_error + header_total

with open('D:\Tools\OpenAeroStruct\openaerostruct\examples\mytest\FunctionalFormatTest_5\Results\ANN_Testing\TestSamplesCombineError.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header_combined)

    # write multiple rows
    writer.writerows(rawDataTestCombine)


AeroPlot, axsA = plt.subplots(1, 2)

# plot the LoD and CL
axsA[0].scatter(outputTestActual[:,0], outputTestPred[:,0])
axsA[0].set_title('LoD_Wing')
axsA[1].scatter(outputTestActual[:,1], outputTestPred[:,1])
axsA[1].set_title('CL_Wing_total')

# plot the Loads
LoadsPlot, axsL = plt.subplots(6, 7)
for i in range(6):
    for j in range(7):
        index = 7 * i + j + 2
        axsL[i, j].scatter(outputTestActual[:,index], outputTestPred[:,index])
        axsL[i, j].set_title(header[index + 32 - 2])

plt.show()