input = "contact_test02.dat";

ContElem =
{
  type = "FiniteStrainContinuum";
  material =
  {
    type = "PlaneStress";
    E    = 1.e6;
    nu   = 0.25;
  };
};

contact =
{
  type = "disc";
  
  radius    = 1.0;
  centre    = [ 8.0 , 1.4 ];
  direction = [0.,-0.5];
  penalty   = 1.0e6;
};

solver =
{
  type     = "NonlinearSolver";
  maxCycle = 40;
  dtime    = 0.2;
};

outputModules = ["vtk","GraphWriter"];

vtk =
{
  type = "MeshWriter";
};

output =
{
  type = "OutputWriter";
  onScreen = true;
};

GraphWriter = 
{
  onScreen = true;
  
  columns = [ "time","disp" ];
  disp = 
  {
    type = "state";
    node = 48;
    dof  = 'v';
  };
};
