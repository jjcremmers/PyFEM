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

input = "two_fibres.dat";

Fibre =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrain";
    E    = 1.e6;
    nu   = 0.25;   
  };
};

Fibre2 =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrain";
    E    = 1.e6;
    nu   = 0.25;   
  };
};

Epoxy =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrain";
    E    = 1.e5;
    nu   = 0.25;   
  };
};

solver =
{
  type     = "LinearSolver";
};

outputModules = ["vtk","graph"];

vtk =
{
  type = "MeshWriter";
};

graph =
{
  type = "GraphWriter";
  onScreen = true;

  columns = [ "disp" , "load" ];

  disp =
  {
    type = "state";
    node = "SinglePoint";
    dof  = "u";
  };

  load =
  {
    type = "fint";
    node = "Right";
    dof  = "u";
  };
};
