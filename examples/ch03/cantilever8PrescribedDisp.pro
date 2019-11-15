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
############################################################################
#  Description: The general PyFEM input file of the example presented in   #
#               section 3.6 of the book, pages 103--105. This input file   #
#               can be used in combination with the standard PyFEM         #
#               executable script PyFEM.py, or in combination with the     #
#               dedicated NewtonRaphson.py script as presented on page 104.#
#                                                                          #
#  Usage:       pyfem cantilever8.pro                                      #
#                 or                                                       #
#               python NewtonRaphson.py cantilever8.pro                    #
############################################################################

input = "cantilever8PrescribedDisp.dat";

ContElem =
{
  type = "FiniteStrainContinuum";

  material =
  {
    type = "PlaneStress";
    E    = 100.0;
    nu   = 0.3;
  };
};

solver =
{
  type = "NonlinearSolver";

  fixedStep = true;
  maxCycle   = 10;
};

outputModules = [ "MeshWriter" , "GraphWriter" ];

GraphWriter = 
{
  onScreen = true;

  columns = [ "disp" , "load" ];

  disp = 
  {
    type = "state";
    node = 48;
    dof  = 'v';
  };
  
  load =
  {
    type = "fint";
    node = 48;
    dof  = 'v';
  };
};
