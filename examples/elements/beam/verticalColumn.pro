input = "verticalColumn.dat";

BeamElem =
{
  type = "BeamNL";
  
  E   = 1.0e5;
  A   = 0.1;
  I   = 1.0e-3;
  G   = 50.0;
  rho = 10.0;
  
  bodyForce = True;
};

solver =
{
  type = "BuckEigSolver";
  
  eigenCount = 5;
};

outputModules = ["vtk","h5"];

vtk =
{
  type = "MeshWriter";

  interval = 1;
};

h5     =
{
  type = "HDF5Writer";
};
