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

RVE =
{
  type = "Prescribed";
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
