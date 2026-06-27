input = "dogbone3D.dat";

ContElem =
{
  type = "PhaseField";

  material =
  {
    type = "Isotropic";
    E    = 2.e6;
    nu   = 0.3;
  };
  
  Gc = 1.0;
  l0 = 1.0;
};

solver =
{
  type = "NonlinearSolver";
  
  solver1 =
  {
    name     = "PhaseField";
    type     = "Linear";
    dofTypes = ["phase"];
  };
  
  solver2 =
  {
    name     = "Displacement";
    type     = "Linear";
    dofTypes = ["u","v","w"];
  };
  
  maxCycle = 50;
};

outputModules = ["vtk","graph"];

vtk =
{
  type = "MeshWriter";
  
  extraFields = "phase";
};

output =
{
  type = "OutputWriter";

  onScreen = true;
};

graph =
{
  type = "GraphWriter";

  columns = [ "disp" , "load" ];

  disp = 
  {
    type = "state";
    node = 10;
    dof  = "u";
  };

  load = 
  {
    type = "fint";
    node = 10;
    dof  = "u";
  };
};
