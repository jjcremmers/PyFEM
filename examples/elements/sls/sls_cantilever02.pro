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

input = "sls_cantilever02.dat";

SLSElem0 =
{
  type = "SLS";

  material = 
  {
    type = "TransverseIsotropic";
    
    E1   = 1.e6;
    E2   = 5.e5;
    nu12 = 0.25;
    G12  = 4.e5;
    rho  = 1.1e3;
  };
  
  theta  = 0.0;     
};

SLSElem1 =
{
  type = "SLS";

  material = 
  {
    type = "TransverseIsotropic";
      
    E1   = 1.e6;
    E2   = 5.e5;
    nu12 = 0.25;
    G12  = 4.e5;
    rho  = 1.2e3;
  };
  
  theta  = 90.0;     
};

SLSElem2 =
{
  type = "SLS";

  material = 
  {
    type = "TransverseIsotropic";
      
    E1   = 1.e6;
    E2   = 5.e5;
    nu12 = 0.25;
    G12  = 4.e5;
    rho  = 1.3e3;
  };
  
  theta  = 0.0;     
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
