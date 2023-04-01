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

from numpy import outer, ones, zeros
from pyfem.materials.MaterialManager import MaterialManager

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
