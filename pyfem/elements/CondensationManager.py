# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers


from numpy import zeros,dot
from scipy.linalg import eigvals,inv

#==============================================================================
#
#==============================================================================

class CondensationManager():

  def __init__( self , condDOF , totDOF ):

    self.condDOF = condDOF
    self.totDOF  = totDOF
    self.intDOF  = self.totDOF - self.condDOF

    self.kwuold    = zeros( shape=( self.condDOF , self.intDOF ) )
    self.kwwinvold = zeros( shape=( self.intDOF  , self.intDOF ) )
    self.fwold     = zeros( self.intDOF  )

    self.kwunew    = zeros( shape=( self.condDOF , self.intDOF ) )
    self.kwwinvnew = zeros( shape=( self.intDOF  , self.intDOF ) )
    self.fwnew     = zeros( self.intDOF )

    self.uold      = zeros( self.condDOF )
    self.wold      = zeros( self.intDOF  )
 
    self.unew      = zeros( self.condDOF )
    self.wnew      = zeros( self.intDOF  )

    self.activeFlag= False
  
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def decondensate( self , elemdat ):

    self.storeU( elemdat.state )

    elemdat.state0 = elemdat.state - elemdat.Dstate
  
    if self.activeFlag:
      tmpArray = dot( self.getKwu().T , elemdat.Dstate )
      tmpArray += self.getFw()

      elemdat.dw = -1.*dot( self.getKwwinv() , tmpArray )
      elemdat.w  = self.getWold()  + elemdat.dw
    else:
      elemdat.dw = zeros( self.intDOF )
      elemdat.w  = zeros( self.intDOF )
    
    self.storeW( elemdat.w )

    elemdat.fullstiff = zeros( shape=( self.totDOF , self.totDOF ) )
    elemdat.fullfint  = zeros( self.totDOF )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def condensate( self , elemdat ):
      
    if self.activeFlag:      
      kwwinv = inv( elemdat.fullstiff[self.condDOF:,self.condDOF:] )
      kwu    = elemdat.fullstiff[self.condDOF:,:self.condDOF]
      fw     = elemdat.fullfint[self.condDOF:]

      self.storeKwwinv( kwwinv )
      self.storeKwu   ( kwu )
      self.storeFw    ( fw )
   
      elemdat.stiff = elemdat.fullstiff[:self.condDOF,:self.condDOF] - \
                        dot( elemdat.fullstiff[:self.condDOF, self.condDOF:] , \
                        dot( kwwinv , kwu ) )
    
      elemdat.fint  = elemdat.fullfint[:self.condDOF] - \
                        dot( elemdat.fullstiff[:self.condDOF,self.condDOF:] , \
                        dot( kwwinv , fw ) )
    else:
      elemdat.stiff = elemdat.fullstiff[:self.condDOF,:self.condDOF]
      elemdat.fint  = elemdat.fullfint[:self.condDOF]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
 
  def commit( self ):

    self.kwuold , self.kwunew = self.kwunew , self.kwuold
    self.kwwinvold , self.kwwinvnew = self.kwwinvnew , self.kwwinvold
    self.fwold , self.fwnew = self.fwnew , self.fwold

    self.uold , self.unew = self.unew , self.uold
    self.wold , self.wnew = self.wnew , self.wold
 
    self.activeFlag = True

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def storeU( self , u ):
    self.unew = u

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def storeW( self , w ):
    self.wnew = w

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def storeKwu( self , kwu ):
    self.kwunew = kwu

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def storeKwwinv( self , kwwinv ):
    self.kwwinvnew = kwwinv

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def storeFw( self , fw ):
    self.fwnew = fw

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getUold( self ):
    return self.uold

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getWold( self ):
    return self.wold

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
  
  def getKwu( self ):
    return self.kwuold

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
    
  def getKwwinv( self ):
    return self.kwwinvold

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getFw( self ):
    return self.fwold
