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
from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from pyfem.util.kinematics      import Kinematics
from numpy import zeros, dot, outer, ones , eye

class SmallStrainContinuum3D( Element ):

  #dofs per element
  dofTypes = [ 'u' , 'v' , 'w' ]
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

  def __type__ ( self ):
    return name

#------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    sData = getElemShapeData( elemdat.coords )

    kin = Kinematics(3,6)
    
    elemdat.outlabel.append(["s11","s22","s33","s23","s13","s12"])

    elemdat.outdata  = zeros( shape=(len(elemdat.nodes),6) )

    for iData in sData:
      
      b = self.getBmatrix( iData.dhdx )

      kin.strain  = dot ( b , elemdat.state )
      kin.dstrain = dot ( b , elemdat.Dstate )
      
      sigma,tang = self.mat.getStress( kin )

      elemdat.stiff += dot ( b.transpose() , dot ( tang , b ) ) * iData.weight
      elemdat.fint  += dot ( b.transpose() , sigma ) * iData.weight

      elemdat.outdata += outer( ones(len(elemdat.nodes)), sigma )
    
    elemdat.outdata *= 1.0 / len(sData) 

  
     
#-------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    kin = Kinematics(3,6)
    
    elemdat.outlabel.append(["s11","s22","s33","s23","s13","s12"])

    elemdat.outdata  = zeros( shape=(len(elemdat.nodes),6) )

    for iData in sData:
      b = self.getBmatrix( iData.dhdx )

      kin.strain  = dot ( b , elemdat.state )
      kin.dstrain = dot ( b , elemdat.Dstate )

      sigma,tang = self.mat.getStress( kin )

      elemdat.fint    += dot ( b.transpose() , sigma ) * iData.weight
      elemdat.outdata += outer( ones(len(self)), sigma )
      
    elemdat.outdata *= 1.0 / len(sData)  

#----------------------------------------------------------------------
    
  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    rho = elemdat.matprops.rho

    for iData in sData:
      N  = self.getNmatrix( iData.h )
      elemdat.mass += dot ( N.transpose() , N ) * rho * iData.weight
     
    elemdat.lumped = sum(elemdat.mass)
   
#--------------------------------------------------------------------------

  def getBmatrix( self , dphi ):

    b = zeros( shape=( 6 , self.dofCount() ) )

    for i,dp in enumerate(dphi):
      b[0,i*3  ] = dp[0]
      b[1,i*3+1] = dp[1]
      b[2,i*3+2] = dp[2]

      b[3,i*3+1] = dp[2]
      b[3,i*3+2] = dp[1]

      b[4,i*3  ] = dp[2]
      b[4,i*3+2] = dp[0]

      b[5,i*3  ] = dp[1]
      b[5,i*3+1] = dp[0]
   
    return b

#------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( 3 , 3*len(h) ) )

    for i,a in enumerate( h ):
      N[0,3*i  ] = a
      N[1,3*i+1] = a
      N[1,3*i+2] = a
   
    return N
