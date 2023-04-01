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
from pyfem.util.kinematics              import Kinematics
from pyfem.elements.SLSgeomdata         import SLSgeomdata
from pyfem.elements.SLSkinematic        import SLSkinematic
from pyfem.elements.SLSutils            import LayerData,SLSparameters,StressContainer
from pyfem.elements.CondensationManager import CondensationManager

from numpy import zeros, dot, outer, ones
 
#==============================================================================
#
#==============================================================================

class SLS( Element ):

  dofTypes = [ 'u' , 'v' , 'w' ]
   
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.param      = SLSparameters( len(elnodes) )
  
    self.kinematic  = SLSkinematic( self.param )
    self.layers     = LayerData( props )

    self.param.nLay = len( self.layers )
    self.stresscont = StressContainer( self.param )

    self.condman    = CondensationManager( self.param.condDOF , self.param.totDOF )
    self.kin        = Kinematics(3,6)
        
  def __type__ ( self ):
    return name
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    self.condman.decondensate( elemdat )

    elemGeomData = SLSgeomdata( elemdat , self.layers )

    self.stresscont.reset()
              
    for iInt,sdat in enumerate(elemGeomData):

      self.kinematic.getDefVecs( sdat , elemdat )

      for iLay,ldat in enumerate(sdat.layerData):

        lamb = ldat.lamb
              
        for iIntZeta,zdat in enumerate(ldat.zetaData):

          bmat = self.kinematic.getBmat( sdat , zdat.zeta , lamb )

          self.kinematic.getStrains( self.kin  , sdat , zdat.zeta , lamb )	 
          self.kin.iMat = ldat.matID

          sigma,tang = self.mat.getStress( self.kin )

          self.stresscont.store( self.mat.outData() , iLay , iIntZeta )          
         
          stiff = dot( bmat.transpose() , dot ( tang , bmat ) )

          self.kinematic.addGeomStiff( stiff , sdat , sigma , lamb , zdat.zeta )
	
          elemdat.fullstiff += stiff * zdat.weight
                    
          elemdat.fullfint  += dot( bmat.transpose() , sigma ) * zdat.weight

    self.appendNodalOutput( self.mat.outLabels() , self.stresscont.getData() )
    
    self.condman.condensate( elemdat )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
   
    self.getTangentStiffness( elemdat )  

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def commit ( self, elemdat ):

    self.condman.commit()
