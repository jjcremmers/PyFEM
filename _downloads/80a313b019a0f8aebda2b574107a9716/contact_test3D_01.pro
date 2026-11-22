input = "contact_test3D_01.dat";

ContElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "Isotropic";
    E    = 1.e6;
    nu   = 0.25;
  };
};

contact =
{
  type = "sphere";
  
  radius    = 1.0;
  centre    = [ 0.0 , 0.0 , 1.0 ];
  direction = [0.,0.,-0.05];
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
