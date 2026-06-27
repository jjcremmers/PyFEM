input = "panel_buckling.dat";

Shell =
{
  type = "ReissnerMindlinShell";

  material =
  {
    E   = 2.1e5;
    nu  = 0.3;
    rho = 7.85e-9;
  };

  thickness = 1.0;

  drillingScale = 1.0e-6;
};

solver =
{
  type = "BuckEigSolver";
};

outputModules = [ "vtk" ];

vtk =
{
  type = "MeshWriter";
  interval = 1;
};