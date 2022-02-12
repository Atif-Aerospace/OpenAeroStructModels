"""
Produce mesh for planform wing
No twist and dihedral
No structural deformation
"""
import numpy as np

def InitialMesh(half_span, kink_location, root_chord, kink_chord, tip_chord, inboard_LE_sweep, outboard_LE_sweep, nx, ny_outboard, ny_inboard):
    # mesh specifications
    # number of chordwise nodal points (should be odd)
    # number of spanwise nodal points for the outboard segment
    # number of spanwise nodal points for the inboard segment
    ny_total = ny_outboard + ny_inboard - 1

    if kink_location == half_span:
        ny_outboard = 1

    # Initialize the 3-D mesh object. Indexing: Chordwise, spanwise, then the 3-D coordinates.
    # We use ny_inboard+ny_outboard-1 because the 2 segments share the nodes where they connect.
    mesh_initial_right = np.zeros((nx, ny_inboard + ny_outboard - 1, 3))


    # For each node,  x, y, and z coordinates.
    # x is streamwise, y is spanwise, and z is up.
    # node for the leading edge at the tip would be specified as mesh[0, 0, :] = np.array([x, y, z]).
    # node at the trailing edge at the root would be mesh[nx-1, ny-1, :] = np.array([x, y, z]).
    # only the right half of the wing is provided due to symmetry.


    ####### THE Z-COORDINATES ######
    # Assume no dihedral, so set the z-coordinate for all the points to 0.
    mesh_initial_right[:, :, 2] = 0.0

    ####### THE Y-COORDINATES ######
    # Using uniform spacing for the spanwise locations of all the nodes within each of the two trapezoidal segments:
    # Outboard
    mesh_initial_right[:, :ny_outboard, 1] = np.linspace(half_span, kink_location, ny_outboard)
    # Inboard
    mesh_initial_right[:, ny_outboard : ny_outboard + ny_inboard, 1] = np.linspace(kink_location, 0, ny_inboard)[1:]

    ###### THE X-COORDINATES ######
    # Start with the leading edge and create some intermediate arrays that we will use
    x_LE = np.zeros(ny_inboard + ny_outboard - 1)

    array_for_inboard_leading_edge_x_coord = np.linspace(0, kink_location, ny_inboard) * np.tan(
        inboard_LE_sweep / 180.0 * np.pi
    )

    array_for_outboard_leading_edge_x_coord = (
        np.linspace(0, half_span - kink_location, ny_outboard) * np.tan(outboard_LE_sweep / 180.0 * np.pi)
        + np.ones(ny_outboard) * array_for_inboard_leading_edge_x_coord[-1]
    )

    x_LE[:ny_inboard] = array_for_inboard_leading_edge_x_coord
    x_LE[ny_inboard : ny_inboard + ny_outboard] = array_for_outboard_leading_edge_x_coord[1:]

    # Then the trailing edge
    x_TE = np.zeros(ny_inboard + ny_outboard - 1)

    array_for_inboard_trailing_edge_x_coord = np.linspace(
        array_for_inboard_leading_edge_x_coord[0] + root_chord,
        array_for_inboard_leading_edge_x_coord[-1] + kink_chord,
        ny_inboard,
    )

    array_for_outboard_trailing_edge_x_coord = np.linspace(
        array_for_outboard_leading_edge_x_coord[0] + kink_chord,
        array_for_outboard_leading_edge_x_coord[-1] + tip_chord,
        ny_outboard,
    )

    x_TE[:ny_inboard] = array_for_inboard_trailing_edge_x_coord
    x_TE[ny_inboard : ny_inboard + ny_outboard] = array_for_outboard_trailing_edge_x_coord[1:]

    # # Quick plot to check leading and trailing edge x-coords
    # plt.plot(x_LE, np.arange(0, ny_inboard+ny_outboard-1), marker='*')
    # plt.plot(x_TE, np.arange(0, ny_inboard+ny_outboard-1), marker='*')
    # plt.show()
    # exit()

    for i in range(0, ny_inboard + ny_outboard - 1):
        mesh_initial_right[:, i, 0] = np.linspace(np.flip(x_LE)[i], np.flip(x_TE)[i], nx)

    # Decompose the matrices into different arrays for processing, printing, plotting, and file writing
    # break def_mesh into x,y,z coordinates
    # mesh_initial_right[0,:,:] leading edge from tip
    # mesh_initial_right[-1,:,:] trailing edge from tip
    # Seems that mesh_undeformed produced by OpenAeroStruct is the left wing while mesh_initial_right is the right wing
    mesh_initial_left = mesh_initial_right.copy()
    mesh_initial_left[:,:,1] = - mesh_initial_right[:,:,1] 
    mesh_initial_left_x = mesh_initial_left[:,:,0]
    mesh_initial_left_y = mesh_initial_left[:,:,1]
    mesh_initial_left_z = mesh_initial_left[:,:,2]

    
    return (ny_total, mesh_initial_right, mesh_initial_left, mesh_initial_left_x, mesh_initial_left_y, mesh_initial_left_z)