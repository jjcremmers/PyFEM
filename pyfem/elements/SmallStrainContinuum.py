################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by                #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since        #
#  then augmented and maintained by Joris J.C. Remmers.                        #
#  All rights reserved.                                                        #
#                                                                              #
#  A github repository, with the most up to date version of the code,          #
#  can be found here:                                                          #
#     https://github.com/jjcremmers/PyFEM/                                     #
#     https://pyfem.readthedocs.io/                                            #	
#                                                                              #
#  The original code can be downloaded from the web-site:                      #
#     http://www.wiley.com/go/deborst                                          #
#                                                                              #
#  The code is open source and intended for educational and scientific         #
#  purposes only. If you use PyFEM in your research, the developers would      #
#  be grateful if you could cite the book.                                     #    
#                                                                              #
#  Disclaimer:                                                                 #
#  The authors reserve all rights but do not guarantee that the code is        #
#  free from errors. Furthermore, the authors shall not be liable in any       #
#  event caused by the use of the program.                                     #
################################################################################

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from pyfem.util.kinematics      import Kinematics
from numpy import zeros, dot, outer, ones , eye

class SmallStrainContinuum( Element ):
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.rank = props.rank

    if self.rank == 2:
      self.dofTypes = [ 'u' , 'v' ]
      self.nstr = 3
      self.outputLabels = ["s11","s22","s12"]
    elif self.rank == 3:
      self.dofTypes = [ 'u' , 'v' , 'w' ]
      self.nstr = 6
      self.outputLabels = ["s11","s22","s33","s23","s13","s12"]

    self.kin = Kinematics(self.rank,self.nstr)

  def __type__ ( self ):
    return name

#------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    sData = getElemShapeData( elemdat.coords )
    
    elemdat.outlabel.append(self.outputLabels)
    elemdat.outdata  = zeros( shape=(len(elemdat.nodes),self.nstr) )

    for iData in sData:
      
      b = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( b , elemdat.state )
      self.kin.dstrain = dot ( b , elemdat.Dstate )
      
      sigma,tang = self.mat.getStress( self.kin )

      elemdat.stiff += dot ( b.transpose() , dot ( tang , b ) ) * iData.weight
      elemdat.fint  += dot ( b.transpose() , sigma ) * iData.weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
     
#-------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    elemdat.outlabel.append(self.outputLabels)
    elemdat.outdata  = zeros( shape=(len(elemdat.nodes),self.nstr) )

    for iData in sData:
      b = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( b , elemdat.state )
      self.kin.dstrain = dot ( b , elemdat.Dstate )

      sigma,tang = self.mat.getStress( self.kin )

      elemdat.fint    += dot ( b.transpose() , sigma ) * iData.weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
 
#-------------------------------------------------------------------------------
    
  def getDissipation ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    for iData in sData:
      b = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( b , elemdat.state )
      self.kin.dstrain = dot ( b , elemdat.Dstate )

      self.mat.getStress( self.kin )

      self.kin.dgdstrain = zeros( 3 )
      self.kin.g = 0.0
            
      elemdat.fint += dot ( b.transpose() , self.kin.dgdstrain ) * iData.weight
      elemdat.diss += self.kin.g * iData.weight
           
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

    b = zeros( shape=( self.nstr , self.dofCount() ) )

    if self.rank == 2:
      for i,dp in enumerate(dphi):
        b[0,i*2+0] = dp[0]
        b[1,i*2+1] = dp[1]
        b[2,i*2+0] = dp[1]
        b[2,i*2+1] = dp[0]
    elif self.rank == 3:
      for i,dp in enumerate(dphi):
        b[0,i*3+0] = dp[0]
        b[1,i*3+1] = dp[1]
        b[2,i*3+2] = dp[2]

        b[3,i*3+1] = dp[2]
        b[3,i*3+2] = dp[1]

        b[4,i*3+0] = dp[2]
        b[4,i*3+2] = dp[0]

        b[5,i*3+0] = dp[1]
        b[5,i*3+1] = dp[0]
   
    return b

#------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( self.rank , self.rank*len(h) ) )

    for i,a in enumerate( h ):
      for j in list(range(self.rank)):
        N[j,self.rank*i+j] = a
    
    return N
