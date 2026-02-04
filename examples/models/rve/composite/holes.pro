input = "holes.dat";

Fibre =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStress";
    E    = 1.e7;
    nu   = 0.25;
  };
};

Epoxy =
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
