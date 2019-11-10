1. The heat transfer simulation is done for an Hemisphere of 144mm diameter with internal heat source generation. The simulation is done for 1 sec with KEpsilon Turbulence model.
2. The solver used is BuoyantBoussinesqPimpleFoam with provided Boundary condition.
3. The folder 0 contains the Boundary conditions of different parameters. We are mainly interested in Temparature change.
4. The folder applications/solvers contains the solver Breastgland (name edited from BuoyantBoussinesqPimpleFoam).
5. The folder constant contains heatSource term (fvOptions), transport Parameters and Turbulence parameters.
6. The folder postProcessing contains the Temparature probe values at different points of the simulations.
7. The folder system contains the mesh, Schemes and Solutions.
