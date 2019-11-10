1. The heat transfer simulation is done for an Hemisphere of 144mm diameter with internal heat source generation. The simulation is done for 1 sec with laminar flow.
2. The solver used is buoyantSimpleFoam with provided Boundary condition.
3. The folder 0 contains the Boundary conditions of different parameters. The BC is applied using GroovyBC library. The parameters of air heat transfer coefficient and ambient temp is added. 
4. The folder constant contains heatSource term (fvOptions), thermophysical properties and Turbulence parameters.
5. The folder system contains the mesh, time parameters, Schemes and Solutions. In time parameter file "controlDict", extra libraries are added which is mandatory for groovyBC along with functions for Temparature probe.

In order to run the simulation, please run blockMesh and checkMesh for checking mesh first. Secondly, type "buoyantSimpleFoam" to run the simulation. 
