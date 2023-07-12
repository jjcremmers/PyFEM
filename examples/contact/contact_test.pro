input = "contact_test.dat";

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
  centre    = [ 5.0 , 2.0 ];
  direction = [0.,-0.5];
  penalty   = 1.0e6;
};

solver =
{
  type     = "NonlinearSolver";
  maxCycle = 20;
  dtime    = 0.1;
};

outputModules = ["vtk"];

vtk =
{
  type = "MeshWriter";
};
