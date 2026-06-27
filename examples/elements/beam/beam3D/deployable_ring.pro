input = "deployable_ring.dat";

# Simplified deployable-ring benchmark inspired by:
# Y. Goto, Y. Watanabe, T. Kasugai, and M. Obata,
# Int. J. Solids Struct. 29(7), 893-909, 1992.
# DOI: 10.1016/0020-7683(92)90024-N

BeamElem =
{
  type = "Beam3D";

  E   = 2.0e5;
  G   = 1.0e5;
  A   = 3.6;
  Ix  = 0.108;
  Iy  = 108.0;
  J   = 108.0;
  rho = 7.85e-9;
  orientation = [ 0. , 0. , 1.0 ];
};

solver =
{
  type = "RiksSolver";

  maxCycle  = 40;
  tol       = 1.0e-3;
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
    node = 24;
    dof  = "rx";
  };

  load =
  {
    type = "fint";
    node = 24;
    dof  = "rx";
  };
};

output =
{
  type = "OutputWriter";
  onScreen = false;
};
