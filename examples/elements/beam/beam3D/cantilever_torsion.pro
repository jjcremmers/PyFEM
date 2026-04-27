input = "cantilever_torsion.dat";

BeamElem =
{
  type = "Beam3D";

  E   = 2.10e5;
  G   = 8.00e4;
  A   = 1.00e-2;
  Ix  = 1.00e-4;
  Iy  = 1.00e-4;
  J   = 2.00e-4;
  rho = 7.85e-9;
  orientation = [ 0.0 , 0.0 , 1.0 ];
};

solver =
{
  type = "NonlinearSolver";

  fixedStep = true;
  maxCycle  = 10;
  tol       = 1.0e-8;
  iterMax   = 20;
};

outputModules = [ "vtk" , "graph" , "output" ];

vtk =
{
  type = "MeshWriter";
  beam = true;
  interval = 1;
};

graph =
{
  type = "GraphWriter";
  onScreen = true;

  columns = [ "disp" , "load" ];

  disp =
  {
    type = "state";
    node = 2;
    dof  = "rx";
  };

  load =
  {
    type = "fint";
    node = 2;
    dof  = "rx";
  };
};

output =
{
  type = "OutputWriter";
  onScreen = true;
};
