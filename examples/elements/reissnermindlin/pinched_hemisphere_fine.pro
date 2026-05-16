############################################################################
#  Reissner-Mindlin shell example: fine pinched hemisphere with 18 degree  #
#  hole, ported from ~/Git/dawn/examples/model/shell/pinched_shell.dat     #
############################################################################

input = "pinched_hemisphere_fine.dat";

Shell =
{
  type = "ReissnerMindlinShell";

  material =
  {
    E   = 6.8253e7;
    nu  = 0.3;
    rho = 1.0;
  };

  thickness = 0.04;

  drillingScale = 1.0e-6;
};

solver =
{
  type = "RiksSolver";

  fixedStep = false;
  maxLam    = 202.0;
  tol       = 1.0e-4;
  iterMax   = 20;
};

outputModules = [ "vtk" , "GraphWriter" ];

vtk =
{
  type = "MeshWriter";
  interval = 1;
};

GraphWriter =
{
  onScreen = true;

  columns = [ "dispx" , "dispy" , "lam" ];

  dispx =
  {
    type = "state";
    node = 1;
    dof  = "u";
  };

  dispy =
  {
    type   = "state";
    node   = 17;
    dof    = "v";
    factor = -1.0;
  };

  lam =
  {
    type = "lam";
  };
};

output =
{
  type = "OutputWriter";
  onScreen = true;
};
