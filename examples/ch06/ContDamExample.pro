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
#  Description: The general PyFEM input file of the shallow truss example  #
#               as shown in the sections 1.6 and 3.2 of the book. In this  #
#               case, the problem is solved using Riks' arclength method.  #
#                                                                          #
#  Usage:       pyfem ContDamExample.pro                                   #
############################################################################

input = "ContDamExample.dat";

ContElem  = 
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrainDamage";

    E    = 100.0;
    nu   = 0.3;
    k    = 1.0;
    kappa0 = 1.0e-6;
    kappac = 1.0e-5;
  };
};

solver =
{
  type = 'NonlinearSolver';

  fixedStep = true;
  maxLam    = 100.0; 
  dtime     = 1.0;
};

outputModules = ["graph","mesh"];

graph =
{
  type = "GraphWriter";

  onScreen = true;

  columns = [ "disp" , "load" ];

  disp =
  {
    type = "state";
    node = 12;
    dof  = 'u';
    factor = 1.0;
  };
  
  load =
  { 
    type = "fint";
    node = 12;
    dof  = 'u';
  };
};

mesh =
{ 
  type = "MeshWriter";
};
