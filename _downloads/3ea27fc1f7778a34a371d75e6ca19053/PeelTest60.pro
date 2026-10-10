input = "PeelTest60.dat";

ContElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrain";
    E    = 100.0;
    nu   = 0.3;
  };
};

InterfaceElem =
{
  type = "Interface";

  material = 
  {
    type = "XuNeedleman";
    
    Tult = 0.5;
    Gc   = 0.1;
  };
};

solver =
{
  type = "DissipatedEnergySolver";

  maxCycle   = 60;
  tol        = 10e-4;
  maxLam     = 50;

  disstype   = "Local";
  switchEnergy = 1.0e-3; 
  maxdTau    = 0.05;
};

outputModules = ["vtk","graph","contour","h5"];

vtk =
{
  type = "MeshWriter";
  
  elementGroup = "ContElem";
};

graph =
{
  type = "GraphWriter";

  onScreen = true;

  columns = ["disp","load"];

  disp =
  {
    type = "state";
    node = 366;
    dof  = "v";
  };
  
  load =
  {
    type = "fint";
    node = 366;
    dof  = "v";
  };
};

contour =
{
  type = "ContourWriter";
  
  nodes = [123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177];
};

h5 =
{
  type = "HDF5Writer";
}
