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
#               section 9.2 of the book, pages 317--321. The amplitude of  #
#               of the sinusoidal imperfection of the beam is h=0.001.     #
#                                                                          #
#  Usage:       pyfem KirchhoffEuler.pro                                   #
############################################################################

input = "KirchhoffEuler.dat";

BeamElem =
{
  type = "KirchhoffBeam";
  E    = 7.2e6;
  A    = 6.0;
  I    = 2.0;
  G    = 5e6;
};

solver =
{
  type = "RiksSolver";

  maxLam    = 35;
  maxFactor = 5;
};

outputModules = [ 'graph' , 'mesh']; 

graph = 
{
  type = "GraphWriter";
  
  onScreen = true;
  
  columns = [ "disp" , "load" ];

  disp =
  { 
    type = "state";
    node = 11;
    dof  = 'v';
  };

  load = 
  {
    type = "fint";
    node = 21;
    dof  = 'u';
    factor = -1.0;
  };
};

mesh =
{
  type = "MeshWriter";
  beam = true;
};

