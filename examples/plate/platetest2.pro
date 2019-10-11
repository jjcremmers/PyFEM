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
#               section 2.6 of the book, pages 53--62.                     #
#                                                                          #
#  Use:         pyfem PatchTest4.pro                                       #
############################################################################

input = "platetest2.dat";

ContElem =
{
  type = "Plate";

  materials = [ "carb" ];

  layers = ["c0"];

  carb =
  {
    type = "PlaneStress";
    E    = 1.e6;
    nu   = 0.25;
    rho  = 1.0e3;
  };

  c0 = 
  {
    material = "carb";
    theta    = 0.;
    thick    = 0.01;
  };

  stack = ["c0"];
};

solver =
{
  type = "LinearSolver";
};

outputModules = ["vtk","output"];

vtk =
{
  type = "MeshWriter";

  interval = 1;
};

output =
{
  type = "OutputWriter";

  onScreen = true;
};
