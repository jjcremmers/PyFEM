input = "beamDyn.dat";

BeamElem =
{
  type = "BeamNL";
  
  E = 1.0e5;
  A = 0.1;
  I = 1.0e-3;
  G = 50.0;
  rho = 10.0;
};

solver =
{
  type = "DynEigSolver";
  
  eigenCount = 5;
};

outputModules = ["vtk"];

vtk =
{
  type = "MeshWriter";

  interval = 1;
};

output =
{
  type = "OutputWriter";

  onScreen = true;
};
