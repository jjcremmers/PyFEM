# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from numpy import outer, ones, zeros
from pyfem.materials.MaterialManager import MaterialManager

class elementData:
  
  outputNames = []
  
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class Element ( list ):

  dofTypes = []

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def __init__ ( self, elnodes , props ):
    list.__init__( self, elnodes )

    self.family     = "CONTINUUM"
    self.history    = {}
    self.current    = {}
    self.solverStat = props.solverStat   

    for name,val in props:
      if name == "material":
        self.matProps = val
        
        self.matProps.rank       = props.rank
        self.matProps.solverStat = self.solverStat
        self.mat = MaterialManager( self.matProps )
      
      setattr( self, name, val )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def dofCount ( self ):

    return len( self ) * len( self.dofTypes )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getNodes ( self ):
    return self

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getType ( self ):
    return self.elemType

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def appendNodalOutput( self , labels , data , weight = 1.0 ):

    for i,name in enumerate(labels):
      if not hasattr( self.globdat , name ):
        self.globdat.outputNames.append( name )

        setattr( self.globdat, name             , zeros( len(self.globdat.nodes) ) )
        setattr( self.globdat, name + 'Weights' , zeros( len(self.globdat.nodes) ) )

      outMat     = getattr( self.globdat , name )
      outWeights = getattr( self.globdat , name + 'Weights' )

      if data.ndim == 1:
        for idx in self.globdat.nodes.getIndices( self ):
          outMat[ idx ]     += data[i]
          outWeights[ idx ] += weight
      else:
        for j,idx in enumerate(self.globdat.nodes.getIndices( self )):
          outMat[ idx ]     += data[j,i]
          outWeights[ idx ] += weight

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def appendElementOutput( self , labels , data , weight = 1.0 ):
      
    if not hasattr( self.globdat , "elementData" ):
      setattr( self.globdat, "elementData" , elementData() )
      
    elemData = getattr( self.globdat , "elementData" )
    
    for i,name in enumerate(labels):
      if not hasattr( elemData , name ):
        elemData.outputNames.append( name )

        setattr( elemData, name             , zeros( len(self.globdat.elements) ) )
        setattr( elemData, name + 'Weights' , zeros( len(self.globdat.elements) ) )

      outMat     = getattr( elemData , name )
      outWeights = getattr( elemData , name + 'Weights' )

      if data.ndim == 1:
        outMat[ self.iElm ]     += data[i]
        outWeights[ self.iElm ] += weight
                   
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def setHistoryParameter ( self, name, val ):
    self.current[name] = val

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getHistoryParameter ( self, name ):
    return self.history[name]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def commitHistory ( self ):
    self.history = self.current.copy()
    self.current = {}

    if hasattr( self , "mat" ):
      self.mat.commitHistory()

  def commit ( self, elemdat ):
    pass
    
#
#
#

  def loadFactor( self ):
    
    return self.solverStat.lam
