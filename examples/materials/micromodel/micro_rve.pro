input = "micro_rve.dat";

MicroElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrainDamage";
    E = 1000.0;
    nu = 0.25;
    k = 1.0;
    kappa0 = 7.5e-4;
    kappac = 2.5e-3;
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
  iterMax = 30;
};
