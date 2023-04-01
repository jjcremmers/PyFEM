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

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros, dot, array, outer
from math import sqrt

class PlaneStrainDamage( BaseMaterial ):

  O1 = array([1.,0.,0.])
  O2 = array([0.,1.,0.])
  O3 = array([0.,0.,1.])

  sc = 1./3.

  def __init__ ( self, props ):

    BaseMaterial.__init__( self, props )

    self.De = zeros( shape = (3,3) )

    self.De[0,0] = self.E*(1.-self.nu)/((1.+self.nu)*(1.-2.*self.nu))
    self.De[0,1] = self.De[0,0]*self.nu/(1.-self.nu)
    self.De[1,0] = self.De[0,1]
    self.De[1,1] = self.De[0,0]
    self.De[2,2] = self.De[0,0]*0.5*(1.-2.*self.nu)/(1.-self.nu)

    self.a1 = (1./(2.*self.k))
    self.a2 = (self.k-1.)/(1.-2.*self.nu)
    self.a3 = 12.*self.k/((1.+self.nu)**2)

    self.a4 = sqrt( self.a2**2+self.a3*self.sc )
    self.O4 = array([self.a4,self.a4,2.*self.a3])

    self.c = self.nu/(self.nu-1.)

    self.setHistoryParameter( 'kappa', 0. )
    self.commitHistory()

    self.outLabels = [ "S11" , "S22" , "S12" , "damage" ]
    self.outData   = zeros(4)

#------------------------------------------------------------------------------
#  pre:  kinematics object containing current strain (kinemtics.strain)
#  post: stress vector and tangent matrix
#------------------------------------------------------------------------------

  def getStress( self, kinematics ):

    kappa = self.getHistoryParameter('kappa')                      

    eps , detadstrain = self.getEquivStrain( kinematics.strain )

    if eps > kappa:
      progDam = True
      kappa   = eps
    else:
      progDam = False
    
    self.setHistoryParameter( 'kappa', kappa )
 
    omega , domegadkappa = self.getDamage( kappa )    

    effStress = dot( self.De , kinematics.strain )

    stress    = ( 1. - omega ) * effStress
    tang      = ( 1. - omega ) * self.De

    if progDam:
      tang += -domegadkappa * outer( effStress , detadstrain )

    self.outData[0:3] = stress
    self.outData[3]   = omega

    return stress , tang           
 
#------------------------------------------------------------------------------
#  pre:  current strain (array of length 3)
#  post: equivalent strain (eps) and its derivative w.r.t. strain array
#------------------------------------------------------------------------------

  def getEquivStrain( self , strain ):                   

    exx = strain[0]
    eyy = strain[1]
    exy = strain[2]
    ezz = self.nu/(self.nu-1.0)*(exx+eyy)

    I1 = exx+eyy+ezz
    J2 = (exx**2+eyy**2+ezz**2-exx*eyy-eyy*ezz-exx*ezz)/3.0+exy**2

    eps = self.a1*(self.a2*I1+sqrt((self.a2*I1)**2+self.a3*J2))

    dexxdstrain = self.O1
    deyydstrain = self.O2
    dexydstrain = 0.5*self.O3
    dezzdstrain = self.c*(dexxdstrain+deyydstrain)

    depsdstrain = zeros(3)

    dI1dstrain = dexxdstrain+deyydstrain+dezzdstrain

    dJ2dstrain  = self.sc*(2.*exx-eyy-ezz)*dexxdstrain
    dJ2dstrain += self.sc*(2.*eyy-exx-ezz)*deyydstrain
    dJ2dstrain += self.sc*(2.*ezz-exx-eyy)*dezzdstrain
    dJ2dstrain += 2.*exy*dexydstrain

    disc = (self.a2*I1)**2+self.a3*J2
    ddiscdI1 = 2.*self.a2**2*I1
    ddiscdJ2 = self.a3

    if disc < 1e-16:
      tmp = 0.
      dtmpdstrain = self.O4

      eps = self.a1*( self.a2*I1 + tmp )

      detadI1 = self.a1*(self.a2)

      detadstrain = detadI1*dI1dstrain + self.a1*dtmpdstrain
    else:   
      tmp = sqrt(disc)
      dtmpdI1 = (0.5/tmp)*ddiscdI1
      dtmpdJ2 = (0.5/tmp)*ddiscdJ2

      eps = self.a1*( self.a2*I1 + tmp )

      detadI1 = self.a1*(self.a2+dtmpdI1)
      detadJ2 = self.a1*(dtmpdJ2)

      detadstrain = detadI1*dI1dstrain + detadJ2*dJ2dstrain

    return eps , depsdstrain

#------------------------------------------------------------------------------
#  pre:  equivalent strain term kappa
#  post: damage (omega) and its derivative w.r.t. kappa (domegadkappa)
#------------------------------------------------------------------------------

  def getDamage( self , kappa ):

    if kappa <= self.kappa0:
      omega        = 0.
      domegadkappa = 0.
    elif self.kappa0 < kappa < self.kappac:
      fac = self.kappac/kappa
      omega        = fac*(kappa-self.kappa0)/(self.kappac-self.kappa0)
      domegadkappa = fac/(self.kappac-self.kappa0)-(omega/kappa)
    else:
      omega        = 1.
      domegadkappa = 0.

    return omega , domegadkappa
