input = "cantilever_plate.dat";

PlateElem =
{
  type = "Plate";

  materials = [ "UD" , "woven" ];

  layers = ["c0","c90","c90","c0"];

  UD =
  {
    E1   = 135000.0;
    E2   = 10000.0;
    nu12 = 0.3;
    G12  = 5000.0;
    rho  = 1600.0;
  };

  woven =
  {
    E1   = 85000.0;
    E2   = 85000.0;
    nu12 = 0.1;
    G12  = 5000.0;
    rho  = 1600.0;
  };

  c0 = 
  {
    material   = "UD";
    theta      = 0.;
    thickness  = 0.25;
  };

  c90 = 
  {
    material   = "UD";
    theta      = 90.;
    thickness  = 0.25;
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
