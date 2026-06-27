input = "macro_fe2.dat";

MacroElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "MicroModel";
    rveFile = "examples/materials/micromodel/micro_rve.pro";
    perturbation = 1.0e-7;
    exportSamples = [0];
  };
};

solver =
{
  type = "NonlinearSolver";
  tol = 1.0e-6;
  iterMax = 12;
  maxCycle = 1;
  loadTable = [1.0];
};

outputModules = ["graph","mesh"];

graph =
{
  type = "GraphWriter";

  columns = ["disp","load"];

  disp =
  {
    type = "state";
    node = 3;
    dof = "u";
  };

  load =
  {
    type = "fint";
    node = 3;
    dof = "u";
  };
};

mesh =
{
  type = "MeshWriter";
};
