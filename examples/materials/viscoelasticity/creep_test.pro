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
# Example: Viscoelastic creep test
#
# This example demonstrates the time-dependent creep behavior of a viscoelastic
# material subjected to constant stress. A single element is loaded with a
# constant force and the displacement response is tracked over time.
#
# The material uses a generalized Maxwell model with multiple relaxation times
# to capture the viscoelastic response.
#

input = "creep_test.dat";

Continuum =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "ViscoElasticity";
    
    # Elastic properties
    E    = 1000.0;         
    nu   = 0.3;           
    Einf = 100.0;         
    
    # Viscoelastic properties - 3 Maxwell elements
    nMaxwell = 3;
    relaxTimes  = [0.1, 1.0, 10.0];   
    relaxModuli = [300.0, 300.0, 300.0];
  };
};

solver =
{
  type = "NonlinearSolver";
  
  maxCycle = 100;
  
  # Time stepping
  tmax = 50.0;
  dtime = 0.5;

  lam = "t*(t<10)+10*(t>=10)";
};

outputModules = ["vtk", "GraphWriter"];

vtk =
{
  type = "MeshWriter";
};

GraphWriter =
{
  onScreen = true;

  columns = ["time", "disp", "stress"];

  time =
  {
    type = "time";
  };

  disp =
  {
    type = "state";
    node = 2;
    dof  = 'u';
  };
 
  stress =
  {
    type = "out";
    node = 1;
    label = "S11";
  };
};
