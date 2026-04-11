input = "clamped_beam3.dat";

ContElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrain";
    E    = 100.0;
    nu   = 0.3;
    rho  = 1.0;
  };
};

solver =
{
  type = "NonlinearSolver";
};

outputModules = ["vtk","rom"];

rom =
{
  type = "ROMSnapshotWriter";
  
  filename = "clamped.h5";
};

vtk =
{
  type = "MeshWriter";
};
