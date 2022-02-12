"""
Produce assumed derformed mesh due to aerodynamics loads
two options
- Control point
- Delta Mesh
"""

import numpy as np
from scipy import interpolate




# Control point base method
# 6 points are used to defined the deformation, 3 on the leading edge while 3 on th trailing edge
# these delta values are normalized with the semi-span
# the first delta is the deformation at root which is always zero
def DeformedMesh_CP(mesh_undeformed_left, fem_origin,
delta_x_LE_1, delta_x_LE_2, delta_x_LE_3,   # Assumed deformation on leading edge: 1, 2, 3 are from tip to root
delta_y_LE_1, delta_y_LE_2, delta_y_LE_3,
delta_z_LE_1, delta_z_LE_2, delta_z_LE_3,
delta_x_TE_1, delta_x_TE_2, delta_x_TE_3,   # Assumed deformation on trailing edge: 1, 2, 3 are from tip to root
delta_y_TE_1, delta_y_TE_2, delta_y_TE_3,
delta_z_TE_1, delta_z_TE_2, delta_z_TE_3 ):
    
    # Leading Edge from ROOT to TIP (Note that the sequence of the number is 3, 2, 1)
    delta_x_LE = np.array([ 0.000000, delta_x_LE_3, delta_x_LE_2, delta_x_LE_1 ])
    delta_y_LE = np.array([ 0.000000, delta_y_LE_3, delta_y_LE_2, delta_y_LE_1 ])
    delta_z_LE = np.array([ 0.000000, delta_z_LE_3, delta_z_LE_2, delta_z_LE_1 ])

    # Trailing Edge from ROOT to TIP (Note that the sequence of the number is 3, 2, 1)
    delta_x_TE = np.array([ 0.000000, delta_x_TE_3, delta_x_TE_2, delta_x_TE_1 ])
    delta_y_TE = np.array([ 0.000000, delta_y_TE_3, delta_y_TE_2, delta_y_TE_1 ])
    delta_z_TE = np.array([ 0.000000, delta_z_TE_3, delta_z_TE_2, delta_z_TE_1 ])

    # number of points for predefined deformation 
    num_delta = delta_x_LE.shape[0]

    # input for spline creation from root to tip
    x_in_spline = (np.linspace(0, 1, num = num_delta))

    # Create spline lines for x, y, and z displacement on Leading edge and trailing edge
    tck_delta_x_LE = interpolate.splrep(x_in_spline, delta_x_LE)
    tck_delta_y_LE = interpolate.splrep(x_in_spline, delta_y_LE)
    tck_delta_z_LE = interpolate.splrep(x_in_spline, delta_z_LE)
    tck_delta_x_TE = interpolate.splrep(x_in_spline, delta_x_TE)
    tck_delta_y_TE = interpolate.splrep(x_in_spline, delta_y_TE)
    tck_delta_z_TE = interpolate.splrep(x_in_spline, delta_z_TE)

    # use the mesh_undeformed_left as a base
    mesh_undeformed_left_x = mesh_undeformed_left[:,:,0]
    mesh_undeformed_left_y = mesh_undeformed_left[:,:,1]
    mesh_undeformed_left_z = mesh_undeformed_left[:,:,2]

    # Displacement at all mesh points
    mesh_delta_left_x = np.zeros(mesh_undeformed_left_x.shape)
    mesh_delta_left_y = np.zeros(mesh_undeformed_left_y.shape)
    mesh_delta_left_z = np.zeros(mesh_undeformed_left_z.shape)

    # Interpolating the deformation along the Leading Edge
    # Note that the mesh_delta_x is defined from tip to root, but no need to use flip as mesh_undeformed_left_y is from tip to root already
    # As the delta_x_LE, delta_y_LE, delta_z_LE, delta_x_TE, delta_y_TE, delta_z_TE, are normalized with the semi-span
    # Here we need to multiply mesh_undeformed_left_y[ 0, 0 ] for leading edge and mesh_undeformed_left_y[-1, 0 ] for trailing edge,
    # as these two values are slightly different from the original semi-span. abs() is used because of the left wing give negative value for y
    mesh_delta_left_x[ 0, : ] = interpolate.splev(mesh_undeformed_left_y[ 0 ]/mesh_undeformed_left_y[ 0, 0 ], tck_delta_x_LE, der = 0) * abs(mesh_undeformed_left_y[ 0, 0 ])
    mesh_delta_left_y[ 0, : ] = interpolate.splev(mesh_undeformed_left_y[ 0 ]/mesh_undeformed_left_y[ 0, 0 ], tck_delta_y_LE, der = 0) * abs(mesh_undeformed_left_y[ 0, 0 ])
    mesh_delta_left_z[ 0, : ] = interpolate.splev(mesh_undeformed_left_y[ 0 ]/mesh_undeformed_left_y[ 0, 0 ], tck_delta_z_LE, der = 0) * abs(mesh_undeformed_left_y[ 0, 0 ])

    # Interpolating the deformation along the Trailing Edge
    mesh_delta_left_x[-1, : ] = interpolate.splev(mesh_undeformed_left_y[-1 ]/mesh_undeformed_left_y[-1, 0 ], tck_delta_x_TE, der = 0) * abs(mesh_undeformed_left_y[-1, 0 ])
    mesh_delta_left_y[-1, : ] = interpolate.splev(mesh_undeformed_left_y[-1 ]/mesh_undeformed_left_y[-1, 0 ], tck_delta_y_TE, der = 0) * abs(mesh_undeformed_left_y[-1, 0 ])
    mesh_delta_left_z[-1, : ] = interpolate.splev(mesh_undeformed_left_y[-1 ]/mesh_undeformed_left_y[-1, 0 ], tck_delta_z_TE, der = 0) * abs(mesh_undeformed_left_y[-1, 0 ])

    # Number of points in the chordwise direction
    nx = mesh_undeformed_left.shape[0]

    # Linear Interpolation of the deformation along each chord using the deformation of the leading adn trailing edges on that chord
    for i in range(1, nx-1):
        aa = mesh_delta_left_x[ 0, : ]
        bb = mesh_delta_left_x[-1, : ]
        mesh_delta_left_x[ i, : ] = aa + i / (nx-1) * (bb - aa)

        aa = mesh_delta_left_y[ 0, : ]
        bb = mesh_delta_left_y[-1, : ]
        mesh_delta_left_y[ i, : ] = aa + i / (nx-1) * (bb - aa)

        aa = mesh_delta_left_z[ 0, : ]
        bb = mesh_delta_left_z[-1, : ]
        mesh_delta_left_z[ i, : ] = aa + i / (nx-1) * (bb - aa)

    # Interpolation
    mesh_delta_left = np.zeros(mesh_undeformed_left.shape)
    mesh_delta_left[:,:,0] = mesh_delta_left_x
    mesh_delta_left[:,:,1] = mesh_delta_left_y
    mesh_delta_left[:,:,2] = mesh_delta_left_z

    # for testing
    f_mesh_delta = 1
    mesh_delta = mesh_delta_left * f_mesh_delta

    # mesh_deformed: mesh after the structural deformation, this mesh will be used in the acutual aerodynamics analysis
    mesh_deformed_left = mesh_undeformed_left + mesh_delta
    mesh_deformed_left_x = mesh_deformed_left[:,:,0]
    mesh_deformed_left_y = mesh_deformed_left[:,:,1]
    mesh_deformed_left_z = mesh_deformed_left[:,:,2]

    # Compute the structural nodes based on the fem_origin location (weighted sum of the LE and TE mesh vertices)
    # Note that here the spar nodes are computed based on the mesh_deformed
    # In aerodynamics no structural analysis are conducted, these points are only used as locations to computeing the loads 
    # s_pts [ny, 3]
    # Taken from load_transfer.py
    # Note that the spar nodes depend ONLY on the LE and TE points
    spar_pts_deformed_left = (1 - fem_origin) * mesh_deformed_left[0, :, :] + fem_origin * mesh_deformed_left[-1, :, :]
    spar_pts_deformed_left_x = spar_pts_deformed_left[:,0]
    spar_pts_deformed_left_y = spar_pts_deformed_left[:,1]
    spar_pts_deformed_left_z = spar_pts_deformed_left[:,2]



    return (mesh_delta_left, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
        mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
        spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z)

# Define the defromed mesh based on mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z
# these will be from the structual design
def DeformedMesh_DM(mesh_undeformed_left, fem_origin, mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z):
    
    
    mesh_delta_left = np.zeros(mesh_undeformed_left.shape)
    mesh_delta_left[:,:,0] = mesh_delta_left_x
    mesh_delta_left[:,:,1] = mesh_delta_left_y
    mesh_delta_left[:,:,2] = mesh_delta_left_z

    # for testing
    f_mesh_delta = 1
    mesh_delta = mesh_delta_left * f_mesh_delta

    # mesh_deformed: mesh after the structural deformation, this mesh will be used in the acutual aerodynamics analysis
    mesh_deformed_left = mesh_undeformed_left + mesh_delta
    mesh_deformed_left_x = mesh_deformed_left[:,:,0]
    mesh_deformed_left_y = mesh_deformed_left[:,:,1]
    mesh_deformed_left_z = mesh_deformed_left[:,:,2]

    # Compute the structural nodes based on the fem_origin location (weighted sum of the LE and TE mesh vertices)
    # Note that here the spar nodes are computed based on the mesh_deformed
    # In aerodynamics no structural analysis are conducted, these points are only used as locations to computeing the loads 
    # s_pts [ny, 3]
    # Taken from load_transfer.py
    # Note that the spar nodes depend ONLY on the LE and TE points
    spar_pts_deformed_left = (1 - fem_origin) * mesh_deformed_left[0, :, :] + fem_origin * mesh_deformed_left[-1, :, :]
    spar_pts_deformed_left_x = spar_pts_deformed_left[:,0]
    spar_pts_deformed_left_y = spar_pts_deformed_left[:,1]
    spar_pts_deformed_left_z = spar_pts_deformed_left[:,2]



    return (mesh_delta_left, # mesh_delta_left_x, mesh_delta_left_y, mesh_delta_left_z,
        mesh_deformed_left, mesh_deformed_left_x, mesh_deformed_left_y, mesh_deformed_left_z,
        spar_pts_deformed_left, spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z)