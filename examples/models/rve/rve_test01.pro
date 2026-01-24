input = "rve_test01.dat";

ContElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStress";
    E    = 1.e6;
    nu   = 0.25;
  };
};

models = ["rve"];

rve =
{
  type = "RVE";
};

solver =
{
  type = "LinearSolver";
};

outputModules = ["vtk"];

vtk =
{
  type = "MeshWriter";
};
