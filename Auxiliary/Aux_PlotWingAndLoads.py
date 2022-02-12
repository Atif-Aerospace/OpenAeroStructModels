"""Plot wing and the Force in the FEM nodes
moments are not shown"""

import matplotlib.pyplot as plt 

def PlotWingAndLoads(mesh_initial_left, mesh_undeformed_left, mesh_deformed_left,
    spar_pts_undeformed_left, spar_pts_deformed_left, loads,
    plotInitial = False, plotUndeformed =True, plotDeformed =True, plotLoads =True):

    nx = mesh_initial_left.shape[0]
    ny_total = mesh_initial_left.shape[1]

    WingMeshPlot = plt.figure(1)
    ax = plt.axes(projection='3d')
    
    if plotInitial == True:
        # plot intial mesh (no twist)
        # spanwise lines from the leading edge
        for i in range(0, nx):
            x_spanwise = mesh_initial_left[i,:,0]
            y_spanwise = mesh_initial_left[i,:,1]
            z_spanwise = mesh_initial_left[i,:,2]
            ax.plot(x_spanwise, y_spanwise, z_spanwise, color='grey')
        # chordwise lines from the tip 
        for i in range(0, ny_total):
            x_chordwise = mesh_initial_left[:,i,0]
            y_chordwise = mesh_initial_left[:,i,1]
            z_chordwise = mesh_initial_left[:,i,2]
            ax.plot(x_chordwise, y_chordwise, z_chordwise, color='grey')

    if plotUndeformed == True:
        # plot mesh after twist and dihedral?, no structural deformation in aerodynamics analysis)
        # spanwise lines from the leading edge
        for i in range(0, nx):
            x_spanwise = mesh_undeformed_left[i,:,0]
            y_spanwise = mesh_undeformed_left[i,:,1]
            z_spanwise = mesh_undeformed_left[i,:,2]
            ax.plot(x_spanwise, y_spanwise, z_spanwise, color='lime')
        # chordwise lines from the tip 
        for i in range(0, ny_total):
            x_chordwise = mesh_undeformed_left[:,i,0]
            y_chordwise = mesh_undeformed_left[:,i,1]
            z_chordwise = mesh_undeformed_left[:,i,2]
            ax.plot(x_chordwise, y_chordwise, z_chordwise, color='lime')
        
        # plot undeformed spar
        spar_pts_undeformed_left_x = spar_pts_undeformed_left[:,0]
        spar_pts_undeformed_left_y = spar_pts_undeformed_left[:,1]
        spar_pts_undeformed_left_z = spar_pts_undeformed_left[:,2]
        ax.plot(spar_pts_undeformed_left_x, spar_pts_undeformed_left_y, spar_pts_undeformed_left_z, color='darkgreen')

    if plotDeformed == True:
        # plot mesh after the structural deformation, this mesh will be used in the acutual aerodynamics analysis
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

        # plot deformed spar
        spar_pts_deformed_left_x = spar_pts_deformed_left[:,0]
        spar_pts_deformed_left_y = spar_pts_deformed_left[:,1]
        spar_pts_deformed_left_z = spar_pts_deformed_left[:,2]
        ax.plot(spar_pts_deformed_left_x, spar_pts_deformed_left_y, spar_pts_deformed_left_z, color='blue')
    
    if plotLoads == True:
        # Plot loads forces
        scale_factor = 50000. # otherwise the arrows will be too long

        loads_Fx = loads[:,0]
        loads_Fy = loads[:,1]
        loads_Fz = loads[:,2]
        loads_Mx = loads[:,3]
        loads_My = loads[:,4]
        loads_Mz = loads[:,5]

        Fx_arrow = loads_Fx / scale_factor
        Fy_arrow = loads_Fy / scale_factor
        Fz_arrow = loads_Fz / scale_factor
        
        for i in range(0, ny_total):
            force_vector_x = [spar_pts_deformed_left_x[i], spar_pts_deformed_left_x[i] + Fx_arrow[i]]
            force_vector_y = [spar_pts_deformed_left_y[i], spar_pts_deformed_left_y[i] + Fy_arrow[i]]
            force_vector_z = [spar_pts_deformed_left_z[i], spar_pts_deformed_left_z[i] + Fz_arrow[i]]
            ax.plot(force_vector_x, force_vector_y, force_vector_z, color='red')

    plt.show()