############################################################################
#  Reissner-Mindlin shell example: curved composite cantilever strip       #
############################################################################

input = "curved_cantilever.dat";

Shell =
{
  type = "ReissnerMindlinShell";

  materials = [ "UD" ];
  layers    = [ "ply0_bot" , "ply90_bot" , "ply90_top" , "ply0_top" ];

  UD =
  {
    E1   = 1.35e5;
    E2   = 1.0e4;
    nu12 = 0.3;
    G12  = 5.0e3;
    G13  = 4.0e3;
    G23  = 3.8e3;
    rho  = 1.6e-9;
  };

  ply0_bot =
  {
    material  = "UD";
    theta     = 0.0;
    thickness = 0.005;
  };

  ply90_bot =
  {
    material  = "UD";
    theta     = 90.0;
    thickness = 0.005;
  };

  ply90_top =
  {
    material  = "UD";
    theta     = 90.0;
    thickness = 0.005;
  };

  ply0_top =
  {
    material  = "UD";
    theta     = 0.0;
    thickness = 0.005;
  };

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
