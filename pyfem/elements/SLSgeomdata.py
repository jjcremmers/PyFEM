############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by            #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since    #
#  then augmented and  maintained by Joris J.C. Remmers.                   #
#  All rights reserved.                                                    #
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

from .Element                   import Element
from pyfem.util.shapeFunctions  import getShapeData
from pyfem.util.kinematics      import Kinematics

from numpy import zeros, dot, sqrt,cos,sin,cross
from numpy import pi
from numpy.linalg import norm,det
from scipy.linalg import eigvals,inv
from scipy.special.orthogonal import p_roots as gauss_scheme

#------------------------------------------------------------------------------
#  Utiliy functions
#------------------------------------------------------------------------------

def unit( a ):
  return a/norm(a)

#------------------------------------------------------------------------------
#  Empty classes
#------------------------------------------------------------------------------

class zetaData:
  pass

#------------------------------------------------------------------------------
#  Empty classes
#------------------------------------------------------------------------------

class layData:
  pass

#------------------------------------------------------------------------------
#  SLSgeomdata
#------------------------------------------------------------------------------

class SLSgeomdata():

  def __init__( self , elemdat , layerData ):

    self.layerData = layerData

    self.initShapeFuncs( elemdat.coords.shape[0] )   
    self.initElement   ( elemdat )

  def __iter__( self ):
    return iter(self.shape)

  def initShapeFuncs( self , nNel ):
    
    '''
    
    '''
    
    if nNel == 8:
      self.shape = getShapeData( 0 , 'Gauss' , 'Quad4' )
    elif nNel == 16:
      self.shape = getShapeData( -1 , 'Gauss' , 'Quad8' )
    else:
      raise NotImplementedError('No SLS element with '+str(nNel)+' nodes available')

    for sdat in self.shape:
      sdat.h    *= 0.5
      sdat.dhdx *= 0.5

    intshape = getShapeData( 0 , 'Gauss' , 'Quad4' )

    for sdat,idat in zip(self.shape,intshape):
      sdat.psi = idat.h
    
    self.zetaSample,self.zetaWeights = gauss_scheme( 2 )

   
  def initElement( self , elemdat ):
  
    '''
    
    '''
  
    e     = zeros( shape=( 3 , 3 ) )
    dd    = zeros( shape=( 3 , 2 ) )
    g0    = zeros( shape=( 3 , 3 ) )
    g1    = zeros( shape=( 3 , 3 ) )    
    econ  = zeros( shape=( 3 , 3 ) )
    lax   = zeros( shape=( 3 , 3 ) ) 
    lamb  = zeros( shape=( 3 , 3 ) ) 
    
    for sdat in self.shape:
    
      height  = -1.0;

      sdat.layerData = []

      for ldat in self.layerData:
        ldatnew          = layData()
        ldatnew.angle    = ldat.angle
        ldatnew.matID    = ldat.matID
        ldatnew.zetaData = []
  
        thick = 2.0 * ldat.thick / self.layerData.totThick
 
        for wght,sample in zip(self.zetaWeights,self.zetaSample):
          zdat = zetaData()
          zdat.isowght = sdat.weight*wght*0.5*thick
          zdat.zeta    = height + 0.5 * thick * ( 1.0 + sample )
          zdat.weight  = 0.
          ldatnew.zetaData.append(zdat)
 
        height += thick 
        sdat.layerData.append(ldatnew)

      sdat.height  = height
    
    midNodes = int(elemdat.coords.shape[0]/2)
    
    cmid   = elemdat.coords[midNodes:,:] + elemdat.coords[:midNodes,:]
    dnodes = elemdat.coords[midNodes:,:] - elemdat.coords[:midNodes,:]
  
    for sdat in self.shape:
      for k in range(2):
        e[:,k]  = dot( sdat.dhdx[:,k] , cmid )
        dd[:,k] = dot( sdat.dhdx[:,k] , dnodes )
   
      e[:,2] = dot( sdat.h , dnodes )

      for i in range(3):
        for j in range(3):
          g0[i,j] = dot( e[:,i] , e[:,j] )
                
      g0inv = inv( g0 )
 
      for i in range(3):
        for j in range(3):           
          econ[i,j] = g0inv[j,0] * e[i,0] + g0inv[j,1] * e[i,1] + g0inv[j,2] * e[i,2]  

      g1[0,0] = 2.0 * dot( e[:,0] , dd[:,0] )
      g1[0,1] = dot( e[:,0] , dd[:,1] ) + dot( e[:,1] , dd[:,0] )
      g1[0,2] = dot(dd[:,0] , e[:,2] )

      g1[1,0] = g1[0,1]
      g1[1,1] = 2.0 * dot( e[:,1] , dd[:,1] )
      g1[1,2] = dot( dd[:,1] ,  e[:,2] )
       
      g1[2,0] = g1[0,2]
      g1[2,1] = g1[1,2]
      g1[2,2] = 0.0

      sdat.gbar = dot ( g0inv , g1 )[:2,:2] 

      T = zeros( shape=(3,3) )
     
      T[:,2] = unit(e[:,2])   
      T[0,0] = 1.0
      T[:,1] = unit(cross(  T[:,2] , T[:,0] ) )   
      T[:,0] = unit(cross(  T[:,1] , T[:,2] ) )

      for ldat in sdat.layerData:
        ldat.lamb = zeros( shape = ( 3 , 3 ) )

        lax[0,0] = cos( ldat.angle )
        lax[1,0] = sin( ldat.angle )

        lax[0,1] = cos( ldat.angle + 0.5*pi )
        lax[1,1] = sin( ldat.angle + 0.5*pi )

        lax[2,2] = 1.0

        lax = dot(T,lax)
         
        for i in range(3):
          for j in range(3):
            ldat.lamb[ i , j ] = sum( econ[:,i] * lax[:,j] )
              
      for ldat in sdat.layerData:
        for zdat in ldat.zetaData:
          gtot = g0 + zdat.zeta*g1
          zdat.weight = zdat.isowght*sqrt(det(gtot))
