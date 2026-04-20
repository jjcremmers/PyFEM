############################################################################
#  Reissner-Mindlin shell example: pinched hemisphere with 18 degree hole  #
############################################################################

input = "pinched_hemisphere.dat";

Shell =
{
  type = "ReissnerMindlinShell";

  material =
  {
    E   = 6.825e7;
    nu  = 0.3;
    rho = 7.85e-9;
  };

  thickness = 0.04;

  drillingScale = 1.0e-6;
};

solver =
{
  type = "RiksSolver";

  fixedStep = false;
  maxCycle  = 20;
  tol       = 1.0e-5;
  iterMax   = 20;
};

outputModules = [ "vtk" , "GraphWriter" , "output" ];

vtk =
{
  type = "MeshWriter";
  interval = 1;
};

GraphWriter =
{
  onScreen = true;

  columns = [ "disp" , "load" ];

  disp =
  {
    type = "state";
    node = 72;
    dof  = "u";
  };

  load =
  {
    type = "fint";
    node = 72;
    dof  = "u";
  };
};

output =
{
  type = "OutputWriter";
  onScreen = true;
};
