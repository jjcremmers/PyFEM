# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from numpy import zeros, dot, outer, ones, eye, sqrt, absolute, linalg,cos,sin,cross
from scipy.linalg import eigvals,inv
from math import pi


#-------------------------------------------------------------------------------
#   class SLSparameters
#-------------------------------------------------------------------------------

class SLSparameters:

  def __init__( self , nNod ):

    if nNod == 8:
      self.totDOF   = 28
      self.condDOF  = 24
      self.nodeDOF  = 3
      self.midNodes = 4
      self.extNodes = 8
      self.intNodes = 4
      self.ansFlag  = True
    elif nNod == 16:
      self.totDOF   = 52
      self.condDOF  = 48
      self.nodeDOF  = 3      
      self.midNodes = 8
      self.extNodes = 16
      self.intNodes = 4      
      self.ansFlag  = False

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def getlam4( lam ):

  '''
  Construct the lam4 operator
  
  Args:
  
  lam  input matrix (3,3)
  
  Output:
  
  lam4
  '''

  lam4 = zeros(shape=(3,3,3,3))

  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          lam4[i,j,k,l]=lam[i,k]*lam[j,l]

  return lam4

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def iso2locbase( iso , lam4 ):

  loc = zeros(6)

  loc[0]=iso[0]*lam4[0,0,0,0]+iso[1]*lam4[1,1,0,0]+iso[2]*lam4[2,2,0,0]+ \
    iso[3]*0.5*(lam4[0,1,0,0]+lam4[1,0,0,0])+ \
    iso[4]*0.5*(lam4[1,2,0,0]+lam4[2,1,0,0])+ \
    iso[5]*0.5*(lam4[2,0,0,0]+lam4[0,2,0,0])
  
  loc[1]=iso[0]*lam4[0,0,1,1]+iso[1]*lam4[1,1,1,1]+iso[2]*lam4[2,2,1,1]+ \
    iso[3]*0.5*(lam4[0,1,1,1]+lam4[1,0,1,1])+ \
    iso[4]*0.5*(lam4[1,2,1,1]+lam4[2,1,1,1])+ \
    iso[5]*0.5*(lam4[2,0,1,1]+lam4[0,2,1,1])

  loc[2]=iso[0]*lam4[0,0,2,2]+iso[1]*lam4[1,1,2,2]+iso[2]*lam4[2,2,2,2]+ \
    iso[3]*0.5*(lam4[0,1,2,2]+lam4[1,0,2,2])+ \
    iso[4]*0.5*(lam4[1,2,2,2]+lam4[2,1,2,2])+ \
    iso[5]*0.5*(lam4[2,0,2,2]+lam4[0,2,2,2])
  
  loc[3]=iso[0]*(lam4[0,0,0,1]+lam4[0,0,1,0])+ \
    iso[1]*(lam4[1,1,0,1]+lam4[1,1,1,0])+ \
    iso[2]*(lam4[2,2,0,1]+lam4[2,2,1,0])+ \
    iso[3]*0.5*(lam4[0,1,0,1]+lam4[0,1,1,0]+ \
		lam4[1,0,0,1]+lam4[1,0,1,0])+ \
    iso[4]*0.5*(lam4[1,2,0,1]+lam4[1,2,1,0]+ \
		lam4[2,1,0,1]+lam4[2,1,1,0])+ \
    iso[5]*0.5*(lam4[2,0,0,1]+lam4[2,0,1,0]+ \
		lam4[0,2,0,1]+lam4[0,2,1,0])
  
  loc[4]=iso[0]*(lam4[0,0,1,2]+lam4[0,0,2,1])+ \
    iso[1]*(lam4[1,1,1,2]+lam4[1,1,2,1])+ \
    iso[2]*(lam4[2,2,1,2]+lam4[2,2,2,1])+ \
    iso[3]*0.5*(lam4[0,1,1,2]+lam4[0,1,2,1]+ \
		lam4[1,0,1,2]+lam4[1,0,2,1])+ \
    iso[4]*0.5*(lam4[1,2,1,2]+lam4[1,2,2,1]+ \
		lam4[2,1,1,2]+lam4[2,1,2,1])+ \
    iso[5]*0.5*(lam4[2,0,1,2]+lam4[2,0,2,1]+ \
		lam4[0,2,1,2]+lam4[0,2,2,1])
  
  loc[5]=iso[0]*(lam4[0,0,2,0]+lam4[0,0,0,2])+ \
    iso[1]*(lam4[1,1,2,0]+lam4[1,1,0,2])+ \
    iso[2]*(lam4[2,2,2,0]+lam4[2,2,0,2])+ \
    iso[3]*0.5*(lam4[0,1,2,0]+lam4[0,1,0,2]+ \
		lam4[1,0,2,0]+lam4[1,0,0,2])+ \
    iso[4]*0.5*(lam4[1,2,2,0]+lam4[1,2,0,2]+ \
		lam4[2,1,2,0]+lam4[2,1,0,2])+ \
    iso[5]*0.5*(lam4[2,0,2,0]+lam4[2,0,0,2]+ \
		lam4[0,2,2,0]+lam4[0,2,0,2])

  return loc

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def iso2loc( iso , lam ):

  if iso.ndim == 1:
    loc = iso2locbase( iso , getlam4( lam ) )
  else:
    loc = iso
    for i,col in enumerate(iso.T):
      loc[:,i] = iso2locbase( col , getlam4( lam ) )

  return loc

#------------------------------------------------------------------------------
#  sigma2oomega
#------------------------------------------------------------------------------

def sigma2omega( sigma , lam ):

  '''
  Transformation of the stresses as obtained from the material manager into the SLS
  local frame of reference.
  
  args:
  
  sigma   Input stress
  lam     deformation operator for the given integration point.
  
  output:
  omega   transformed stress
  '''

  omega = zeros(6)

  omega[0]=lam[0,0]*lam[0,0]*sigma[0]+ \
    lam[0,1]*lam[0,1]*sigma[1]+ \
    lam[0,2]*lam[0,2]*sigma[2]+ \
    2*(lam[0,0]*lam[0,1]*sigma[3])+ \
    2*(lam[0,1]*lam[0,2]*sigma[4])+ \
    2*(lam[0,0]*lam[0,2]*sigma[5])
      
  omega[1]=lam[1,0]*lam[1,0]*sigma[0]+ \
    lam[1,1]*lam[1,1]*sigma[1]+ \
    lam[1,2]*lam[1,2]*sigma[2]+ \
    2*(lam[1,0]*lam[1,1]*sigma[3])+ \
    2*(lam[1,1]*lam[1,2]*sigma[4])+ \
    2*(lam[1,0]*lam[1,2]*sigma[5])
  
  omega[2]=lam[2,0]*lam[2,0]*sigma[0]+ \
    lam[2,1]*lam[2,1]*sigma[1]+ \
    lam[2,2]*lam[2,2]*sigma[2]+ \
    2*(lam[2,0]*lam[2,1]*sigma[3])+ \
    2*(lam[2,1]*lam[2,2]*sigma[4])+ \
    2*(lam[2,0]*lam[2,2]*sigma[5])

  omega[3]=lam[0,0]*lam[1,0]*sigma[0]+ \
    lam[0,0]*lam[1,1]*sigma[3]+ \
    lam[0,0]*lam[1,2]*sigma[5]+ \
    lam[0,1]*lam[1,0]*sigma[3]+ \
    lam[0,1]*lam[1,1]*sigma[1]+ \
    lam[0,1]*lam[1,2]*sigma[4]+ \
    lam[0,2]*lam[1,0]*sigma[5]+ \
    lam[0,2]*lam[1,1]*sigma[4]+ \
    lam[0,2]*lam[1,2]*sigma[2]
    
  omega[4]=lam[1,0]*lam[2,0]*sigma[0]+ \
    lam[1,0]*lam[2,1]*sigma[3]+ \
    lam[1,0]*lam[2,2]*sigma[5]+ \
    lam[1,1]*lam[2,0]*sigma[3]+ \
    lam[1,1]*lam[2,1]*sigma[1]+ \
    lam[1,1]*lam[2,2]*sigma[4]+ \
    lam[1,2]*lam[2,0]*sigma[5]+ \
    lam[1,2]*lam[2,1]*sigma[4]+ \
    lam[1,2]*lam[2,2]*sigma[2]

  omega[5]=lam[0,0]*lam[2,0]*sigma[0]+ \
    lam[0,0]*lam[2,1]*sigma[3]+ \
    lam[0,0]*lam[2,2]*sigma[5]+ \
    lam[0,1]*lam[2,0]*sigma[3]+ \
    lam[0,1]*lam[2,1]*sigma[1]+ \
    lam[0,1]*lam[2,2]*sigma[4]+ \
    lam[0,2]*lam[2,0]*sigma[5]+ \
    lam[0,2]*lam[2,1]*sigma[4]+ \
    lam[0,2]*lam[2,2]*sigma[2]

  return omega

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class Layer: 
  pass

#------------------------------------------------------------------------------
#  LayerData
#------------------------------------------------------------------------------

class LayerData:

  def __init__( self , props ):

    self.layers   = []
    self.totThick = 0.

    if hasattr( props , "layers" ):
      for layID in props.layers:
        layprops = getattr( props , layID )
       
        layer          =  Layer()
        layer.thick    =  layprops.thickness
        layer.theta    =  layprops.theta*pi/180
        
        if props.material.type == "MultiMaterial":      
          layer.matID  =  props.material.materials.index(layprops.material)          
          matprops = getattr(props.material,layprops.material)
          layer.rho    =  matprops.rho
        else:
          layer.matID  = 0
          layer.rho    =  props.material.rho
                              
        self.totThick += layprops.thickness

        self.layers.append( layer )
    else:
      layer         = Layer()
      layer.thick   = 1.0
      
      if hasattr( props , "theta" ):
        layer.theta = props.theta
      else:
        layer.theta = 0.0
        
      layer.matID   = 0
      
      if hasattr( props.material , "rho" ):
        layer.rho = props.material.rho
        
      self.totThick = 1.0
      self.layers.append( layer )
      
  def __iter__( self ):
    return iter( self.layers )  

  def __len__( self ):
    return len(self.layers)


#-------------------------------------------------------------------------------
#  StressContainer 
#-------------------------------------------------------------------------------


class StressContainer:

  def __init__( self , param ):
    
    '''
    
    '''
    
    self.nLay = param.nLay
    self.nMid = param.midNodes
    self.nNod = param.extNodes
    self.reset()

  def reset( self , nStr = 6 ):
  
    '''
    Erases all outputdata and sets weights to zero.
    '''
    
    self.nStr    = nStr
      
    self.data    = zeros( shape = ( self.nLay , self.nStr , self.nNod ) )
    self.weights = zeros( self.nLay )
    

  def store( self , matData , iLay , iIntZeta ):
  
    '''
    Stores the material data in the correct array. In case of a 1 layer element, the
    material data in the bottom integration points is stored in the bottom nodes and
    the data for the top integration points in the top nodes. In case of a 
    multi-layered model, all material data is stored with separate labels for all nodes.
    '''
          
    if self.nStr != len(matData):
      self.reset( len(matData) )
      
    if self.nLay == 1:
      if iIntZeta == 0:
        self.data[ 0,:,:self.nMid] += outer( matData , ones(self.nMid) )       
        self.weights[ 0 ] += 1
      elif iIntZeta == 1:
        self.data[ 0 , : , self.nMid: ] += outer( matData , ones(self.nMid) )            
    else:
      self.data[ iLay , : , : ] += outer( matData , ones(self.nNod) )
      self.weights[ iLay ] += 1

  def getData( self ):
  
    '''
    Returns the weight averaged data of all nodes.
    '''

    for iLay in range(self.nLay):
      self.data[iLay,:,:] *= 1.0/self.weights[iLay]
    return self.data.reshape(self.nLay*self.nStr,self.nNod).T

  def getLabels( self ):
  
    '''
    Returns the labels of the data inthe container.
    '''
      
    origlabel = ["s11","s22","s33","s13","s23","s12"]
    if self.nLay == 1:
      return origlabel
    else:
      labels = []
      for iLay in range(self.nLay):
        for ll in origlabel:
          labels.append( "lay"+str(iLay)+"-"+ll)
      return labels
