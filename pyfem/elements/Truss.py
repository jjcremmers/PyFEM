# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element
from pyfem.util.transformations import toElementCoordinates, toGlobalCoordinates

from numpy import zeros, dot, array, eye, outer
from scipy.linalg import norm

class Truss ( Element ):

  #Number of dofs per element
  dofTypes = ['u','v']
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    #Initialize the history parameter
    self.setHistoryParameter( 'sigma', 0. )
    self.commitHistory()
    
    self.family = "BEAM"

  def getTangentStiffness ( self, elemdat ):

    a  = toElementCoordinates( elemdat.state , elemdat.coords )
    Da = toElementCoordinates( elemdat.Dstate, elemdat.coords )
    a0 = a - Da

    #Compute l0 and strains

    self.l0 = norm( elemdat.coords[1]-elemdat.coords[0] )
   
    ( epsilon , Depsilon ) = self.getStrain( a , a0 )

    #Compute the stress increment (multiplied
    #with the undeformed cross-sectional area)
    Dsigma = elemdat.props.E * Depsilon

    #Compute the current stress (multiplied
    #with the undeformed cross-sectional area)
    sigma = self.getHistoryParameter('sigma') + Dsigma

    #Update the history parameter
    self.setHistoryParameter( 'sigma', sigma )

    #Compute BL in the element coordinate system

    BL = self.getBL( a )

    #Compute the element stiffness in the element coordinate system
    KL = elemdat.props.E * elemdat.props.Area * self.l0 * outer( BL , BL )

    KNL = self.getKNL( sigma , elemdat.props.Area )

    elStiff = KL + KNL

    #Rotate element tangent stiffness to the global coordinate system
    elemdat.stiff = toGlobalCoordinates( elStiff , elemdat.coords )

    #Compute the element internal force vector in the element coordinate system
    elFint = self.l0 * sigma * elemdat.props.Area * BL

    #Rotate element fint to the global coordinate system
    elemdat.fint = toGlobalCoordinates( elFint , elemdat.coords )

#-----------------------------------------------------------------

  def getInternalForce ( self, elemdat ):

    #Compute the current state vector

    a  = toElementCoordinates( elemdat.state , elemdat.coords )
    Da = toElementCoordinates( elemdat.Dstate, elemdat.coords )

    a0 = a - Da

    #Compute l0 and strains
    self.l0 = norm( elemdat.coords[1]-elemdat.coords[0] )

    ( epsilon , Depsilon ) = self.getStrain( a , a0 )
  
    #Compute the stress increment (multiplied
    #with the undeformed cross-sectional area)
    Dsigma = elemdat.props.E * Depsilon

    #Compute the current stress (multiplied
    #with the undeformed cross-sectional area)
    sigma = self.getHistoryParameter('sigma') + Dsigma

    #Update the history parameter
    self.setHistoryParameter( 'sigma', sigma )

    #Compute BL in the parent element coordinate system
    BL = self.getBL( a )

    #Compute the element internal force vector in the element coordinate system
    elFint = self.l0 * sigma * elemdat.props.Area * BL

    #Rotate element fint to the global coordinate system
    elemdat.fint = toGlobalCoordinates( elFint, elemdat.coords )

#------------------------------------------

  def getStrain( self , a , a0 ):

    epsilon  = (a[2]-a[0])/self.l0 + 0.5*((a[2]-a[0])/self.l0)**2 + 0.5*((a[3]-a[1])/self.l0)**2
    epsilon0 = (a0[2]-a0[0])/self.l0 + 0.5*((a0[2]-a0[0])/self.l0)**2 + 0.5*((a0[3]-a0[1])/self.l0)**2

    #Compute the strain increment
    Depsilon = epsilon -epsilon0

    return epsilon,Depsilon

#-------------------------------------------

  def getBL( self , a ):

    BL = zeros( 4 )

    BL[0] = (-1./self.l0)*(1.+(a[2]-a[0])/self.l0)
    BL[1] = (-1./self.l0)*(a[3]-a[1])/self.l0
    BL[2] = -BL[0]
    BL[3] = -BL[1]

    return BL
#-------------------------------------------

  def getKNL( self , sigma , A0 ):

    KNL = zeros( (4,4) )

    KNL[:2,:2] =  (sigma * A0/self.l0)*eye(2)
    KNL[:2,2:] = -(sigma * A0/self.l0)*eye(2)
    KNL[2:,:2] = KNL[:2,2:]
    KNL[2:,2:] = KNL[:2,:2]

    return KNL
