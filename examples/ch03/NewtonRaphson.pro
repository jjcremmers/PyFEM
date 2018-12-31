
############################################################################
#  Description: Simulation of a clamped beam with a point load             #
#                                                                          #
#  Usage:       pyfem NewtonRaphson.pro                                    #
############################################################################

input = "NewtonRaphson.dat";

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

solver =
{
  type = "NonlinearSolver";

  maxCycle = 20;
};

outputModules = ["mesh" , "graph" ]; 

mesh =
{
  type     = "MeshWriter";
  prefix   = "NewtonRaphson";
  interval = 5;
};

graph =
{
  type = "GraphWriter";
  onScreen = true;

  columns = [ "disp" , "load" ];

  disp = 
  {
    type = "state";
    node = 21;
    dof  = "v";
  };

  load = 
  {
    type = "lam";
  };
};
