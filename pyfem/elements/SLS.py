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

from .Element import Element
from pyfem.util.kinematics       import Kinematics
from pyfem.elements.SLSgeomdata  import SLSgeomdata
from pyfem.elements.SLSkinematic import SLSkinematic
from pyfem.elements.SLSutils     import LayerData,SLSparameters,StressContainer
from pyfem.elements.CondensationManager     import CondensationManager

from numpy import zeros, dot, outer, ones
 
#==============================================================================
#
#==============================================================================

class SLS( Element ):

  dofTypes = [ 'u' , 'v' , 'w' ]
   
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    param = SLSparameters( len(elnodes) )
  
    self.kinematic  = SLSkinematic( param )
    self.layers     = LayerData( props )

    param.nLay      = len( self.layers )
    self.stresscont = StressContainer( param )

    self.condman    = CondensationManager( 24 ,28 )
    self.kin        = Kinematics(3,6)

        
  def __type__ ( self ):
    return name
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    self.condman.decondensate( elemdat )

    self.slsgeomdata = SLSgeomdata( elemdat , self.layers )

    elemdat.outlabel.append( self.stresscont.getLabels() )

    self.stresscont.reset()

    wei = 0.0
      
    for iInt,sdat in enumerate(self.slsgeomdata):

      self.kinematic.getDefVecs( sdat , elemdat )

      for iLay,ldat in enumerate(sdat.layerData):

        lamb = ldat.lamb
              
        for iIntZeta,zdat in enumerate(ldat.zetaData):

          bmat = self.kinematic.getBmat( sdat , zdat.zeta , lamb )

          self.kinematic.getStrains( self.kin  , sdat , zdat.zeta , lamb )	 
          self.kin.iMat = ldat.matID

          # calculate stresses and tangent	      
          sigma,tang = self.mat.getStress( self.kin )

          # calulate linear part of stiffness matrix        
          tmpMat = dot( bmat.T , dot ( tang , bmat ) )

          # add non-linear part stiffness matrix
          self.kinematic.addGeomStiff( tmpMat , sdat , sigma , lamb , zdat.zeta )
	
          elemdat.fullstiff += tmpMat * zdat.weight
          
          #construct internal forces vector
          elemdat.fullfint  += dot( bmat.transpose() , sigma ) * zdat.weight

          #store stresses
          #elemdat.outdata += outer( ones(len(elemdat.nodes)), sigma*0.125 )

          self.stresscont.store( sigma , iLay , iIntZeta )

          wei += zdat.weight

    elemdat.outdata = self.stresscont.getStress()

#    print(elemdat.outdata)

    self.condman.condensate( elemdat )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
   
    self.getTangentStiffness( elemdat )  

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    rho = elemdat.matprops.rho

    for iData in sData:
      N  = self.getNmatrix( iData.h )
      elemdat.mass += dot ( N.transpose() , N ) * rho * iData.weight
     
    elemdat.lumped = sum(elemdat.mass)

  def commit ( self, elemdat ):

    self.condman.commit()
