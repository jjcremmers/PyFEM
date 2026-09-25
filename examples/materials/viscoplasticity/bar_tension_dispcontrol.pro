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

input = "bar_tension_dispcontrol.dat";

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
    syield = 10.0;       
    hard   = 0.0;      
    
    # Viscoplastic properties
    gamma  = 0.1;        
    n      = 5.0;          
  };
};

solver =
{
  type = "NonlinearSolver";
  iterMax = 50;

  dtime = 0.0005;

  loadTable = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0];
};

outputModules = ["vtk", "GraphWriter"];

vtk =
{
  type = "MeshWriter";
};

GraphWriter =
{
  onScreen = true;

  columns = ["disp", "load"];

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
