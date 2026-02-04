input = "cantilever_plate.dat";

PlateElem =
{
  type = "Plate";

  materials = [ "aluminium" ];

  layers = ["c0"];

  aluminium =
  {
    E    = 7.0e4;
    nu   = 0.3;
    rho  = 2400.0;
  };

  c0 = 
  {
    material   = "aluminium";
    theta      = 0.;
    thickness  = 1.0;
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
