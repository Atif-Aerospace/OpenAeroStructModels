"""produce file for plot in three.js"""
import numpy as np
import sys
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(threshold=sys.maxsize)

def ThreeJSFile(fileName, mesh_deformed_left, mesh_point_force, mesh_delta_left):
    
    nx = mesh_point_force.shape[0] - 1  # from LE to TE
    ny = mesh_point_force.shape[1] - 1   # from Tip to Root

    wingMesh = []

    FX = []
    FY = []
    FZ = []

    deltaX = []
    deltaY = []
    deltaZ = []

    
    # get the locations of 4 points for each panel
    for i in range(0,nx):
        for j in range(0,ny):
            # Point 1 (i,j)
            p1 = np.copy(mesh_deformed_left[i,j,:])
            # Point 2 (i+1,j)
            p2 = np.copy(mesh_deformed_left[i+1,j,:])
            # Point 3 (i+1,j+1)
            p3 = np.copy(mesh_deformed_left[i+1,j+1,:])
            # Point 4 (i,j+1)
            p4 = np.copy(mesh_deformed_left[i,j+1,:])

            # two triangular surface
            triangular1 = np.concatenate((p1, p2, p3))
            triangular2 = np.concatenate((p1, p3, p4))

            # add the wingMesh to array
            wingMesh = np.append(wingMesh, triangular1)
            wingMesh = np.append(wingMesh, triangular2)

            # add the panel FX to array
            FX = np.append(FX, mesh_point_force[i,j,0])
            FX = np.append(FX, mesh_point_force[i+1,j,0]) 
            FX = np.append(FX, mesh_point_force[i+1,j+1,0])
            
            FX = np.append(FX, mesh_point_force[i,j,0])
            FX = np.append(FX, mesh_point_force[i+1,j+1,0])
            FX = np.append(FX, mesh_point_force[i,j+1,0]) 

            # add the panel FY to array
            FY = np.append(FY, mesh_point_force[i,j,1])
            FY = np.append(FY, mesh_point_force[i+1,j,1]) 
            FY = np.append(FY, mesh_point_force[i+1,j+1,1])
            
            FY = np.append(FY, mesh_point_force[i,j,1])
            FY = np.append(FY, mesh_point_force[i+1,j+1,1])
            FY = np.append(FY, mesh_point_force[i,j+1,1]) 

            # add the panel FZ to array
            FZ = np.append(FZ, mesh_point_force[i,j,2])
            FZ = np.append(FZ, mesh_point_force[i+1,j,2]) 
            FZ = np.append(FZ, mesh_point_force[i+1,j+1,2])
            
            FZ = np.append(FZ, mesh_point_force[i,j,2])
            FZ = np.append(FZ, mesh_point_force[i+1,j+1,2])
            FZ = np.append(FZ, mesh_point_force[i,j+1,2]) 


            # add the panel deltaX to array
            deltaX = np.append(deltaX, mesh_delta_left[i,j,0])
            deltaX = np.append(deltaX, mesh_delta_left[i+1,j,0]) 
            deltaX = np.append(deltaX, mesh_delta_left[i+1,j+1,0])
            
            deltaX = np.append(deltaX, mesh_delta_left[i,j,0])
            deltaX = np.append(deltaX, mesh_delta_left[i+1,j+1,0])
            deltaX = np.append(deltaX, mesh_delta_left[i,j+1,0]) 

            # add the panel deltaY to array
            deltaY = np.append(deltaY, mesh_delta_left[i,j,1])
            deltaY = np.append(deltaY, mesh_delta_left[i+1,j,1]) 
            deltaY = np.append(deltaY, mesh_delta_left[i+1,j+1,1])
            
            deltaY = np.append(deltaY, mesh_delta_left[i,j,1])
            deltaY = np.append(deltaY, mesh_delta_left[i+1,j+1,1])
            deltaY = np.append(deltaY, mesh_delta_left[i,j+1,1]) 

            # add the panel deltaZ to array
            deltaZ = np.append(deltaZ, mesh_delta_left[i,j,2])
            deltaZ = np.append(deltaZ, mesh_delta_left[i+1,j,2]) 
            deltaZ = np.append(deltaZ, mesh_delta_left[i+1,j+1,2])
            
            deltaZ = np.append(deltaZ, mesh_delta_left[i,j,2])
            deltaZ = np.append(deltaZ, mesh_delta_left[i+1,j+1,2])
            deltaZ = np.append(deltaZ, mesh_delta_left[i,j+1,2]) 


    templateName = "ThreeJsPlot\models\AeroTemplate.json"
    templateFile  = open(templateName, "r")
    lines = templateFile.readlines()
    templateFile.close()

    #  write wing mesh to file
    lines[12] = "            \"array\": " + np.array2string(wingMesh, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x}) + '\n'

    #  write FX to file
    lines[18] = "            \"array\": " + np.array2string(FX, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x}) + '\n'

    #  write FY to file
    lines[24] = "            \"array\": " + np.array2string(FY, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x}) + '\n'

    #  write FZ to file
    lines[30] = "            \"array\": " + np.array2string(FZ, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x}) + '\n'
    
    #  write deltaX to file
    lines[36] = "            \"array\": " + np.array2string(deltaX, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x}) + '\n'

    #  write deltaY to file
    lines[42] = "            \"array\": " + np.array2string(deltaY, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x}) + '\n'

    #  write deltaZ to file
    lines[48] = "            \"array\": " + np.array2string(deltaZ, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x}) + '\n'
    
    writeFile  = open(fileName, "w")
    writeFile.writelines(lines)
    writeFile.close()

    # FZString = np.array2string(FZ, separator=', ', formatter={'float_kind':lambda x: "%.6f" % x})
    # f.write(FZString + "\n")