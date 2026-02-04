input = "twofibres.dat";

Fibre =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStress";
    E    = 1.e5;
    nu   = 0.25;
  };
};

Epoxy =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStress";
    E    = 1.e4;
    nu   = 0.3;
  };
};

models = ["rve"];

rve =
{
  type = "RVE";
  unitStrain = [0.034,0.15,0.01];
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
