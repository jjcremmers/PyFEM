input = "clamped_beam.dat";

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
  type = "DynEigSolver";
};

outputModules = ["vtk"];

vtk =
{
  type = "MeshWriter";
};
