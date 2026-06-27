input = "micro_rve_fibres.dat";

Fibre =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStress";
    E = 1.0e5;
    nu = 0.22;
  };
};

Epoxy =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStress";
    E = 1.0e4;
    nu = 0.30;
  };
};

models = ["rve"];

rve =
{
  type = "RVE";
  boundaryType = "Periodic";
  unitStrain = [0.0,0.0,0.0];
};

solver =
{
  type = "NonlinearSolver";
  tol = 1.0e-8;
  iterMax = 25;
};
