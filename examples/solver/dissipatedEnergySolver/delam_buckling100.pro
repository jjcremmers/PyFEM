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
#  Description: The delamination buckling example (section 3.1) from       #
#               the journal paper (100 elements in x-direction):           #
#                                                                          #
#               E.Borjesson, J.J.C. Remmers, M. Fagerstrom (2022)          #
#               A generalised path-following solver for robust analysis    #
#               of material failure, Computational Mechanics               #
#               doi: 10.1007/s00466-022-02175-w                            #
#                                                                          #
#  Usage:       pyfem delam_buckling200.pro                                #
############################################################################

input = "delam_buckling100.dat";

ContElem =
{
  type = "FiniteStrainContinuum";

  material =
  {
    type = "PlaneStrain";
    E    = 135e3;
    nu   = 0.18;
  };
};

InterfaceElem =
{
  type = "Interface";

  material = 
  {
    type = "XuNeedleman";
    
    Tult = 75.0;
    Gc   = 0.2;
  };
};

solver =
{
  type = "DissipatedEnergySolver";

  maxCycle   = 30;
  tol        = 10e-4;
  maxLam     = 5000;
  lam        = 1.0;

  disstype   = "Local";
  switchEnergy = 1.0e-3; 
  maxdTau    = 1.0;
};

outputModules = ["vtk","graph"];

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
    type   = "state";
    node   = 202;
    dof    = "v";
    factor = -1.0;
  };
    
  load =
  {
    type   = "fint";
    node   = 202;
    dof    = "u";
    factor = -1.0;
  };
};
