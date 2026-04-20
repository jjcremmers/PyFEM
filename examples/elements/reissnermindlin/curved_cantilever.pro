############################################################################
#  Reissner-Mindlin shell example: curved cantilever strip                 #
############################################################################

input = "curved_cantilever.dat";

Shell =
{
  type = "ReissnerMindlinShell";

  material =
  {
    E   = 2.1e5;
    nu  = 0.3;
    rho = 7.85e-9;
  };

  thickness = 0.02;

  drillingScale = 1.0e-6;
};

solver =
{
  type = "NonlinearSolver";

  fixedStep = true;
  maxCycle  = 20;
  tol       = 1.0e-6;
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
    node = 5;
    dof  = "w";
  };

  load =
  {
    type = "fint";
    node = 5;
    dof  = "w";
  };
};

output =
{
  type = "OutputWriter";
  onScreen = true;
};
