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
  boundaryType = "Prescribed";
  unitStrain = [0.034,0.15,0.01];
};

solver =
{
  type = "LinearSolver";
};

outputModules = ["vtk","graph"];

vtk =
{
  type = "MeshWriter";
};

graph =

{
  type = "GraphWriter";


}
