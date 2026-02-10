input = "thermal_test01.dat";

ContElem =
{
  type = "ThermoContinuum";

  material =
  {
    heatConductivity = 100.0;
    heatCapacity = 100.0;
  };
};

BoundElem =
{
  type = "ThermalBC";
  
  h = 5.0;
  extTemp = 25.0;
};

solver =
{
  type = "NonlinearSolver";
  loadFunc = t*(t<10)+10*(t>=10);
  maxCycle = 40;
};

outputModules = ["vtk"];

vtk =
{
  type = "MeshWriter";
};

