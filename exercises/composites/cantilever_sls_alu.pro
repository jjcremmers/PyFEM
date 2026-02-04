
input = "cantilever_sls.dat";

SLSElem =
{
  type = "SLS";

  material = 
  {
    type = "Isotropic";
    E    = 7.e4;
    nu   = 0.3;
    rho  = 2.4e3;
  };
};

solver =
{
  type = "LinearSolver";
};

outputModules = ["vtk","output"];

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
