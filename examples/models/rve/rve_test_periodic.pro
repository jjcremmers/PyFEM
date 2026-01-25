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
  boundaryType = "Periodic";
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

  columns = [ "sxx" , "exx" ];
  
  sxx = 
  {
    type = "equivStress";
    comp = 0;
  };
  
  exx =
  {
    type = "equivStrain";
    comp = 0;
  };
};
