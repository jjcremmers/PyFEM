# Tutorial: The dissipated energy pathfollowing solver

This example is taken from the paper:

E. Borjesson, J.J.C. Remmers, M. Fagerstrom (2022)  
A generalised path-following solver for robust analysis of material failure, Computational Mechanics  
https://doi.org/10.1007/s00466-022-02175-w

In section 3.1 a delamination buckling simulation is considered. The input file for this case can be found in the directory `examples/solver/dissipatedEnergSolver`.

- `delam_buckling100.pro`: simulation with 450 elements (400 finite strain continuum elements and 50 interface elements)
- `delam_buckling200.pro`: simulation with 800 continuum and 50 interface elements

These input files can be executed by typing:
```bash
pyfem delam_buckling200.pro
```

The result load–displacement curve shows snap-back behavior typical of delamination buckling. Run the example and inspect the solver output or plot files written to the working directory.
