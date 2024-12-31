################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
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
from pyfem.materials.MatUtils     import vonMisesStress,hydrostaticStress,Hardening
from numpy import zeros, ones, dot, array, outer
from math import sqrt

class IsotropicHardeningPlasticity( BaseMaterial ):

  def __init__ ( self, props ):

    self.tolerance = 1.0e-6

    BaseMaterial.__init__( self, props )

    self.syield0 = self.syield
    self.ebulk3  = self.E / ( 1.0 - 2.0*self.nu )
    self.eg2     = self.E / ( 1.0 + self.nu )
    self.eg      = 0.5*self.eg2
    self.eg3     = 3.0*self.eg
    self.elam    = ( self.ebulk3 - self.eg2 ) / 3.0

    self.hardLaw = Hardening( props )

    self.ctang = zeros(shape=(6,6))

    self.ctang[:3,:3] = self.elam

    self.ctang[0,0] += self.eg2
    self.ctang[1,1] = self.ctang[0,0]
    self.ctang[2,2] = self.ctang[0,0]

    self.ctang[3,3] = self.eg
    self.ctang[4,4] = self.ctang[3,3]
    self.ctang[5,5] = self.ctang[3,3]
 
    self.setHistoryParameter( 'sigma'  , zeros(6) )
    self.setHistoryParameter( 'eelas'  , zeros(6) )
    self.setHistoryParameter( 'eplas'  , zeros(6) )
    self.setHistoryParameter( 'eqplas' , zeros(1) )

    self.commitHistory()

    #Set the labels for the output data in this material model
    self.outLabels = [ "S11" , "S22" , "S33" , "S23" , "S13" , "S12" , "Epl" ]
    self.outData   = zeros(7)

#------------------------------------------------------------------------------
#  pre:  kinematics object containing current strain (kinemtics.strain)
#  post: stress vector and tangent matrix
#------------------------------------------------------------------------------

  def getStress( self, kinematics ):

    eelas  = self.getHistoryParameter('eelas')   
    eplas  = self.getHistoryParameter('eplas')   
    eqplas = self.getHistoryParameter('eqplas')   
    sigma  = self.getHistoryParameter('sigma') 
  
    tang = self.ctang

    eelas += kinematics.dstrain

    sigma += dot( self.ctang , kinematics.dstrain )

    smises = vonMisesStress( sigma )

    syield , hard = self.hardLaw.getHardening( eqplas )

    print(syield,eqplas)

    if smises > ( 1.0 + self.tolerance ) * syield:

      shydro = hydrostaticStress( sigma )
   
      flow = sigma

      flow[:3] = flow[:3]-shydro*ones(3)
      flow *= 1.0/smises

      syield = self.syield0

      deqpl = 0.0
      rhs   = syield

      k = 0
    
      while abs(rhs) > self.tolerance * self.syield0:

        k = k+1

        if k > 10:
          printf("rrrr")

        rhs   = smises-self.eg3*deqpl - syield
        deqpl = deqpl+rhs/(self.eg3+hard)

        syield , hard = self.hardLaw.getHardening( eqplas + deqpl )
 
        print("SS",syield,eqplas+deqpl)

      eplas[:3] +=  1.5 * flow[:3] * deqpl
      eelas[:3] += -1.5 * flow[:3] * deqpl

      eplas[3:] +=  3.0 * flow[:3] * deqpl
      eelas[3:] += -3.0 * flow[:3] * deqpl

      sigma = flow * syield
      sigma[:3] += shydro * ones(3)

      eqplas += deqpl

      #output = 0.5*deqpl*(self.syield0+syield)
     
      effg   = self.eg*syield / smises
      effg2  = 2.0*effg
      effg3  = 3.0*effg
      efflam = 1.0/3.0 * ( self.ebulk3-effg2 )
      effhdr = self.eg3 * self.hard/(self.eg3+self.hard)-effg3
     
      tang[:3,:3] = efflam
    
      for i in range(3):
        tang[i,i]     += effg2
        tang[i+3,i+3] += effg

      tang += effhdr*outer(flow,flow)
 
    self.setHistoryParameter( 'eelas' , eelas  )
    self.setHistoryParameter( 'eplas' , eplas  )
    self.setHistoryParameter( 'sigma' , sigma  )
    self.setHistoryParameter( 'eqplas', eqplas )

    # Store output eplas

    self.outData[:6] = sigma
    self.outData[6]  = eqplas

    return sigma , tang           
 
