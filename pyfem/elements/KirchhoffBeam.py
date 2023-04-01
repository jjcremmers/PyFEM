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

from .Element import Element
from pyfem.util.transformations import getRotationMatrix

from numpy import zeros, dot, array, eye, outer, mat, empty,sqrt
from scipy.linalg import norm
from math import atan2, sin, cos, tan

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class KirchhoffBeam ( Element ):

  #dofs per element
  dofTypes = [ 'u' , 'v' , 'rz' ]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.EA = self.E * self.A
    self.EI = self.E * self.I

    self.intpoints = zeros(3)
    self.weights   = zeros(3)
    
    self.weights[0] = 5./9.
    self.weights[1] = 8./9.
    self.weights[2] = 5./9.
        
    self.intpoints[0] = -sqrt(3./5.)
    self.intpoints[1] =  0.
    self.intpoints[2] =  sqrt(3./5.)
    
    self.family = "BEAM"

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def __type__ ( self ):
    return name

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    l0  = norm( elemdat.coords[2]-elemdat.coords[0] )
    jac = 0.5 * l0
    
    a_bar = self.glob2Elem( elemdat.state , elemdat.coords )
            
    fint  = zeros(9);
    stiff = zeros( elemdat.stiff.shape ); 
           
    for xi,alpha in zip( self.intpoints , self.weights ):
      
      bu = self.getBu( l0 , xi )
      bw = self.getBw( l0 , xi )
      c  = self.getC ( l0 , xi )
      
      epsl = dot( bu , a_bar ) + 0.5*(dot( bw , a_bar ) )**2
      chi  = dot( c  , a_bar )
             
      N = self.EA * epsl
      M = self.EI * chi
              
      wght = jac * alpha
    
      fint  += N * bu * wght
      fint  += ( N * dot( bw , a_bar ) * bw + M * c ) * wght
    
      stiff += self.EA * outer( bu , bu ) * wght
      stiff += self.EA * dot( bw , a_bar ) * outer( bu , bw ) * wght
      stiff += self.EA * dot( bw , a_bar ) * outer( bw , bu ) * wght
      stiff += ( self.EI * outer( c , c ) + \
                 self.EA * (dot( bw , a_bar ))**2 * outer( bw , bw ) + \
                 N  * outer( bw , bw ) ) * wght
                         
    stiff[4,4] = 1.0
    stiff[5,5] = 1.0                     
    
    elemdat.fint  = self.elem2Glob( fint  , elemdat.coords )
    elemdat.stiff = self.elem2Glob( stiff , elemdat.coords )    

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
    
  def getInternalForce ( self, elemdat ):

    EA = elemdat.props.E * elemdat.props.A
    EI = elemdat.props.E * elemdat.props.I
    l0 = norm( elemdat.coords[2]-elemdat.coords[0] )

    intpoints = zeros(3)
    weights   = zeros(3)
    
    weights[0] = 5./9.
    weights[1] = 8./9.
    weights[2] = 5./9.
        
    intpoints[0] = -sqrt(3./5.)
    intpoints[1] =  0.
    intpoints[2] =  sqrt(3./5.)
    
    a_bar = elemdat.state
            
    for xi,alpha in zip(intpoints,weights):
      
      bu = self.getBu( l0 , xi )
      bw = self.getBw( l0 , xi )
      c  = self.getC ( l0 , xi )
      
      epsl = dot( bu , a_bar ) + 0.5*(dot( bw , a_bar ) )**2
      chi  = dot( c  , a_bar )
             
      N    = EA * epsl
      M    = EI * chi
              
      wght = 0.5 * l0 * alpha
    
      elemdat.fint  += N * bu * wght
      elemdat.fint  += ( N * dot( bw , a_bar ) * bw + M * c ) * wght
  
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBu( self , l0 , xi ):

    Bu = zeros( 9 )

    Bu[0] = -1.0/l0
    Bu[3] = -4.0*xi/l0
    Bu[6] =  1.0/l0
  
    return Bu

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBw( self , l0 , xi ):

    Bw = zeros(9)

    Bw[1] = 1.5/l0*(xi*xi-1.0)
    Bw[2] = 0.25*(3*xi*xi-2.0*xi-1.0)
    Bw[7] = -1.5/l0*(xi*xi-1.0)
    Bw[8] = 0.25*(3*xi*xi+2.0*xi-1.0)
        
    return Bw

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getC( self , l0 , xi ):

    C = zeros( 9 )

    l02 = l0*l0
    
    C[1] = 6.0*xi/l02
    C[2] = (3.0*xi-1.0)/l0
    C[7] = -6.0*xi/l02
    C[8] = (3.0*xi+1.0)/l0
       
    return C

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def glob2Elem( self , a , coords ):

    a_bar = empty( a.shape )

    R     = eye( 9)
    crd   = zeros( shape=(2,2) )
        
    crd[0,:] = coords[0,:]
    crd[1,:] = coords[2,:]

    R[:2,:2]   = getRotationMatrix( crd )
    R[6:8,6:8] = R[:2,:2]

    a_bar      = dot( R , a )

    return a_bar

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def elem2Glob( self , a_bar , coords ):

    a = empty( a_bar.shape )
    R     = eye( 9)
    crd   = zeros( shape=(2,2) )
    
    crd[0,:] = coords[0,:]
    crd[1,:] = coords[2,:]
    
    R[:2,:2]   = getRotationMatrix( crd )
    R[6:8,6:8] = R[:2,:2]
    
    if len(a_bar.shape) == 1:
      a        = dot( R.transpose() , a_bar )
    elif len(a_bar.shape) == 2:
      a        = dot( R.transpose(), dot ( a_bar , R ) )
    return a

