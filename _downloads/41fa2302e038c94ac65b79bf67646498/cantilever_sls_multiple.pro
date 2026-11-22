

input = "cantilever_sls_multiple.dat";

SLSElem0 =
{
  type = "SLS";

  material = 
  {
    type = "TransverseIsotropic";
    
    E1   = 135000.0;
    E2   = 10000.0;
    nu12 = 0.3;
    G12  = 5000.0;
    rho  = 1600.0;
  };
  
  theta  = 0.0;     
};

SLSElem90 =
{
  type = "SLS";

  material = 
  {
    type = "TransverseIsotropic";
    
    E1   = 135000.0;
    E2   = 10000.0;
    nu12 = 0.3;
    G12  = 5000.0;
    rho  = 1600.0;
  };
  
  theta  = 90.0;     
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
