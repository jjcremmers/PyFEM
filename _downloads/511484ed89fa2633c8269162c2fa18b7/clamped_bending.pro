input = "clamped_bending.dat";

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
  type = "LinearSolver";
};

outputModules = [ "vtk" , "output" ];

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