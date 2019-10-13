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
from pyfem.util.dataStructures import Properties
from pyfem.util.dataStructures import GlobalData

from numpy import zeros, array
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness


from numpy import zeros, array
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness

#################################
# Step 1:                       #
# Initialize (Delta a)          #
#################################




class DissipationNRGSolver:

  def __init__( self , props , globdat ):

    self.tol = 1.0e-4

    print "\n  Starting dissipation solver ..\n"

  def run( self , props , globdat ):

    globdat.cycle += 1
    
    dofCount = len(globdat.dofs)

    a  = globdat.state
    Da = globdat.Dstate

    Da[:] = zeros( dofCount )
    fint  = zeros( dofCount ) 
    fext  = zeros( dofCount ) 

    print '================================='
    print ' Load step %i' % globdat.cycle
    print '================================='
    print '  NR iter : L2-norm residual'
     
    #fext = fext + Dfext
  
    globdat.iiter = 0 

    K = assembleTangentStiffness( props, globdat )
    
    error = 1.
    
    while error > self.tol:

      globdat.iiter += 1
	      
      da = globdat.dofs.solve( K, fext-fint )

      Da[:] += da[:]
      a [:] += da[:]

      fint = assembleInternalForce   ( props, globdat )
      K    = assembleTangentStiffness( props, globdat )
    
      error = globdat.dofs.norm( fext-fint )               
    
      print '  Iter', globdat.iiter, ':', error

      if globdat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')

    # Converged
    
    elements.commitHistory()

    if globdat.cycle == 10:
      globdat.active = False 







      # Solve for new displacement vector, load factor      

      if self.method == 'force-controlled':
        da = dofs.solve( K, lam*fhat-fint )
   
      elif self.method == 'nrg':
        h  =  0.5 * lam0 * fhat
        w  = -0.5 * dot ( (a-Da) , fhat )
        g  =  0.5 * dot ( ( lam0 * Da - Dlam * ( a[:] - Da[:] ) ) , fhat ) - dtau
  
        d1 = dofs.solve( K , lam*fhat - fint )
        d2 = dofs.solve( K , -1.0*fhat )

        denom  = dot ( h , d2 ) - w

        da     = d1 - ( d2 * ( dot( h , d1 ) + g ) ) / denom
        dlam   = -g - ( dot( -1.0*h , d1 ) - g * ( 1.0 + denom ) ) / denom;
 
        Dlam   += dlam
        lam    += dlam
      else:
        raise RuntimeError('Method not known')
   
      # Update displacements

      Da[:] += da[:]
      a [:] += da[:]

      # Solve for new displacement vector, load factor      
  
      K    = assembleTangentStiffness( props, globdat )
      fint = assembleInternalForce( props, globdat )
    
      # Check convergence

      error  = globdat.dofs.norm( lam*fhat-fint ) / globdat.dofs.norm( lam*fhat )

      # Increment the Newton-Raphson iteration counter
      # and print error

      globdat.iiter += 1

      print '   Iter %5i  : %4.2e ' %(iiter,error)

      if globdat.iiter == iterMax:
        plotCurve( output )  
        raise RuntimeError('Newton-Raphson iterations did not converge!')

    # If converged, calculate the amount of energy that has been dissipated in the \
    # previous step.

    print '---------------------------------'
    print '   C O N V E R G E D'
  
    if method == 'force-controlled':
      print '   Diss. nrg   : %1.3e ' %dissnrg

    Da[:]  = zeros( len(dofs) )
    lam0 = lam
  
    if method == 'force-controlled':
      if dissnrg > switchnrg:
        print '   Switch to nrg diss. arc-length'
        self.method = 'nrg'
        Dlam   = 0.;
        dtau   = 0.25*switchnrg
      else:
        lam   += Dlam
    else:
      Dlam = 0.;
      dtau *= pow(0.5,0.25*(iiter-5))
      if dtau > maxnrg:
        dtau = maxnrg

    globdat.elements.commitHistory()

    globdat.dissnrg = 0.5 * dot( ( lam0 * Da - (lam-lam0) * ( a - Da ) ) , fhat )
 
