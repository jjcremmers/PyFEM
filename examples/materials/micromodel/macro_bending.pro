input = "macro_bending.dat";

MacroElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "MicroModel";
    rveFile = "micro_rve_fibres.pro";
    perturbation = 1.0e-7;
  };
};

solver =
{
  type = "NonlinearSolver";
  tol = 1.0e-6;
  iterMax = 15;
  maxCycle = 1;
  loadTable = [1.0];
};

outputModules = ["graph","mesh"];

graph =
{
  type = "GraphWriter";

  columns = ["utop","ubot","ftop","fbot"];

  utop =
  {
    type = "state";
    node = 8;
    dof = "u";
  };

  ubot =
  {
    type = "state";
    node = 2;
    dof = "u";
  };

  ftop =
  {
    type = "fint";
    node = 8;
    dof = "u";
  };

  fbot =
  {
    type = "fint";
    node = 2;
    dof = "u";
  };
};

mesh =
{
  type = "MeshWriter";
};
