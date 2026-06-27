############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stable version can be downloaded from the web-site:          #
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/jjcremmers/PyFEM                                  #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #  
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################

#
# Example: Rate-dependent viscoplastic bar in tension
#
# This example demonstrates the rate-dependent behavior of a viscoplastic
# material. A notched bar is subjected to tension with different loading rates.
# The viscoplastic model shows that higher loading rates result in higher
# peak stresses before plastic flow occurs.
#
# The Perzyna overstress model is used to capture the rate dependency.
#

input = "bar_tension.dat";

Continuum =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "ViscoPlasticity";
    
    # Elastic properties
    E      = 200000.0;     
    nu     = 0.3;         
    
    # Yield properties
    syield = 250.0;       
    hard   = 1000.0;      
    
    # Viscoplastic properties
    gamma  = 0.001;        
    n      = 1.0;          
  };
};

solver =
{
  type = "NonlinearSolver";
  
  maxCycle = 50;
  
  # Time stepping for rate-dependent loading
  tmax = 1.0;
  dtime = 0.02;
};

outputModules = ["vtk", "GraphWriter"];

vtk =
{
  type = "MeshWriter";
};

GraphWriter =
{
  onScreen = true;

  columns = ["time", "disp", "load", "plastic"];

  time =
  {
    type = "time";
  };

  disp =
  {
    type = "state";
    node = 17;
    dof  = 'u';
  };
 
  load =
  {
    type = "fint";
    node = load_nodes;
    dof  = 'u';
  };
  
  plastic =
  {
    type = "out";
    node = 8;
    label = "EqPl";
  };
};
