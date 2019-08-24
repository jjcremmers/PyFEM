############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stabke version can be downloaded from the web-site:          #
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
from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array, pi
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness, assembleMassMatrix

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class DynEigSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol     = 1.0e-3
    self.iterMax = 10 

    BaseModule.__init__( self , props )

    self.fext  = zeros( len(globdat.dofs) )  
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):
      
    K,fint  = assembleTangentStiffness( props, globdat )
         
    M,mlump = assembleMassMatrix      ( props , globdat )

    eigenvals , eigenvecs = globdat.dofs.eigensolve( K , M )

    globdat.state = eigenvecs
  
    globdat.elements.commitHistory()

    globdat.active = False 

    self.printResults( eigenvals )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printResults( self , eigenvals):

    print('\n======================================')
    print(' eigenfrequencies')
    print('======================================')
    print(' Mode  Eigen   Freq')
    
    for i,f in enumerate(eigenvals):
      print(' %3i : %6.4e rad/s  %6.4e Hz' %(i+1,f,f/(2.0*pi)))
      
    print('======================================\n')
