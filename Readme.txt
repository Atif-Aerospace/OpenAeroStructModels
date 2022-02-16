---------------------------------- v4 -----------------------------------
change to Wimpress parameterization
change wing_area to wing_area_geo
change input/output files for addtional variables

---------------------------------- v5 -----------------------------------
start to implement sampling and surrogate
aero surrogate trained with tensorflow
initial implementation of structures in the functional format


---------------------------------- v6 -----------------------------------
clean unused files
more in and out for the structure function
Aero validation done: same with the provided example


---------------------------------- v7 -----------------------------------
Combined Main for iterative process
    - SingleRunCombinedIterative.py

Implementation of Neural Networks with Tensoflow
    - TrainSurrogateAero_Test_1.py
    - folder NeuralNetworks

aero
- add drag polar
    - Aerodynamics\Aero_AnalysisComponent.py
    
- add one more option for using delta mesh as inputs to aero model
    - Aerodynamics\Aero_IntegratedWorkflow.py
    - Aerodynamics\Aero_DerformedMesh.py
    - SingleRunAero_DM.py

- add mesh point force as output of aero for plotting the lift contour

- Loads doubled at root: in v6 this is performed in Aux_PlotWingAndLoads now it's moved to Aero_AnalysisComponent
    - Aerodynamics\Aero_AnalysisComponent.py
    - Auxiliary\Aux_PlotWingAndLoads.py

struct
- Loads become a matrix input
    - Structures\Struct_IntegratedWorkflow.py
    - Structures\Struct_InputFileWimpress.py
    - SingleRunStructure.py

Add three.js for plotting


fix mistake 
# mesh_delta_left = mesh_deformed_left - wingInfo["mesh"] # wingInfo["mesh"] is the initial mesh
mesh_delta_left = mesh_deformed_left - mesh_undeformed_left
- Structures\Struct_AnalysisComponent.py

fix bug arrays values should be defined using np.copy rather just "=", otherise the values will be changed unexpectedly
e.g. in the drag polar
    - Aerodynamics\Aero_AnalysisComponent.py
    - Structures\Struct_AnalysisComponent.py

Note: In SingleRunCombinedIterative.py
we need to change both aif and sif to make the palnform value the same
sif.wing_area_wpref, sif.aspect_ratio_wpref may not be equal to aif.wing_area_wpref, aif.aspect_ratio_wpref

---------------------------------- v8 -----------------------------------

move drag polar out of the aerodynamic analysis function: this is because in iterative approach, the drag polar should be produced after convergence
    - Aerodynamics\Aero_AnalysisComponent.py	
    - Aerodynamics\Aero_IntegratedWorkflow.py
    - SingleRunAero_DM.py

