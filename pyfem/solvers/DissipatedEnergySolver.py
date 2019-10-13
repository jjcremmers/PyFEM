############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
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
from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array, dot
from pyfem.fem.Assembly import assembleTangentStiffness

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class DissipatedEnergySolver( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol       = 1.0e-4
    self.optiter   = 5
    self.iterMax   = 10
    self.maxdTau   = 1.0e20

    self.factor    = 1.0
    self.maxLam    = 1.0e20

    dofCount    = len(globdat.dofs)

    BaseModule.__init__( self , props )

    self.method    = "force-controlled"
    self.Dlam      = 1.0

    globdat.lam    = 1.0
    globdat.dTau   = 0.0

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def run( self , props , globdat ):

    globdat.cycle += 1
   
    a    = globdat.state
    Da   = globdat.Dstate
    fhat = globdat.fhat
 
    self.printHeader( globdat.cycle )
      
    error         = 1.
    globdat.iiter = 0
    lam0          = globdat.lam

    K,fint = assembleTangentStiffness( props, globdat )  

    while error > self.tol:

      if self.method == 'force-controlled':
        da = globdat.dofs.solve( K, globdat.lam*fhat-fint )

      elif self.method == 'nrg-controlled':
        h  =  0.5 * lam0 * fhat
        w  = -0.5 * dot ( (a-Da) , fhat )
        g  =  0.5 * dot ( ( lam0 * Da - self.Dlam * ( a[:] - Da[:] ) ) , fhat ) - globdat.dtau
  
        d1 = globdat.dofs.solve( K , globdat.lam*fhat - fint )
        d2 = globdat.dofs.solve( K , -1.0*fhat )

        denom  = dot ( h , d2 ) - w

        da     = d1 - ( d2 * ( dot( h , d1 ) + g ) ) / denom
        dlam   = -g - ( dot( -1.0*h , d1 ) - g * ( 1.0 + denom ) ) / denom;
 
        self.Dlam   += dlam
        globdat.lam += dlam
     
      else:
        raise RuntimeError('Method not known')
   
      # Update displacements

      Da[:] += da[:]
      a [:] += da[:]

      # Solve for new displacement vector, load factor      
  
      K,fint = assembleTangentStiffness( props, globdat )
    
      # Check convergence

      error  = globdat.dofs.norm( globdat.lam*fhat-fint ) / globdat.dofs.norm( globdat.lam*fhat )

      # Increment the Newton-Raphson iteration counter
      # and print error

      globdat.iiter += 1

      self.printIteration( globdat.iiter , error )

    # If converged, calculate the amount of energy that has been dissipated in the \
    # previous step.

    dissnrg = 0.5 * dot( ( lam0 * Da - (globdat.lam-lam0) * ( a - Da ) ),fhat )

    self.printConverged( globdat.iiter , dissnrg )

    Da[:]  = zeros( len(globdat.dofs) )
  
    if self.method == 'force-controlled':
      if dissnrg > self.switchEnergy:
        print('   Switch to nrg diss. arc-length')
        self.method       = 'nrg-controlled'
        globdat.dtau = 0.25*self.switchEnergy
      else:
        globdat.lam += self.Dlam
    else:
      self.Dlam = 0.
      globdat.dtau *= pow(0.5,0.25*(globdat.iiter-self.optiter))
      if globdat.dtau > self.maxdTau:
        globdat.dtau = self.maxdTau
    
    globdat.elements.commitHistory()

    globdat.fint = fint

    if globdat.lam > self.maxLam or globdat.cycle > self.maxCycle:
      globdat.active=False


#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printHeader( self , cycle):

    print('\n======================================')
    print(' Load step %i' % cycle)
    print('======================================')
    print('  iter # : L2-norm residual')

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printIteration( self , iiter , error ):

    print('   %5i : %4.2e ' %(iiter,error))

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printConverged( self , iiter , dissnrg ):

    print('--------------------------------------')
    print(' Converged in %i iterations' %iiter)
    if self.method == 'force-controlled':
      print(' Dissipated energy : %1.3e ' %dissnrg)
