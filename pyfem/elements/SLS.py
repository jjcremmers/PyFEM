# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element
from pyfem.util.kinematics              import Kinematics
from pyfem.elements.SLSgeomdata         import SLSgeomdata
from pyfem.elements.SLSkinematic        import SLSkinematic
from pyfem.elements.SLSutils            import LayerData,SLSparameters,StressContainer
from pyfem.elements.CondensationManager import CondensationManager
from pyfem.util.shapeFunctions  import getElemShapeData

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
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getMassMatrix ( self, elemdat ):

    elemGeomData = SLSgeomdata( elemdat , self.layers )
             
    for sdat in elemGeomData:
      for ldat in sdat.layerData:              
        rho = ldat.rho
        for zdat in ldat.zetaData:
        
          H = self.kinematic.getHmat( sdat , zdat.zeta )            
          elemdat.mass += dot ( H.T, H ) * rho * zdat.weight
     
    elemdat.lumped = sum(elemdat.mass)
   
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( self.rank , self.rank*len(h) ) )

    for i,a in enumerate( h ):
      for j in list(range(self.rank)):
        N[j,self.rank*i+j] = a
    
    return N  
    
