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

class TimoshenkoBeam ( Element ):

  #dofs per element
  dofTypes = [ 'u' , 'v' , 'rz' ]

  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )
    self.family = "BEAM"

  def __type__ ( self ):
    return name

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    EA = elemdat.props.E * elemdat.props.A
    EI = elemdat.props.E * elemdat.props.I
    GA = 5./6.*elemdat.props.G * elemdat.props.A
    
    l0 = norm( elemdat.coords[2]-elemdat.coords[0] )
    
    intpoints = zeros(2)
    weights   = zeros(2)
    
    weights[0] = 1.
    weights[1] = 1.
       
    intpoints[0] = -0.577
    intpoints[1] =  0.577
    
    a_bar = self.toElementCoordinates( elemdat.state , elemdat.coords )
                
    fint  = zeros(9);
    stiff = zeros( elemdat.stiff.shape ); 
           
    for xi,alpha in zip(intpoints,weights):
      
      ht = self.getHt( xi )
      bu = self.getBu( l0 , xi )
      bw = self.getBw( l0 , xi )
      bt = self.getBt( l0 , xi )
      
      eps  = dot( bu , a_bar ) + 0.5*(dot( bw , a_bar ) )**2
      gam  = dot( ht , a_bar ) + dot( bw , a_bar )
      chi  = dot( bt , a_bar )
             
      N    = EA * eps
      Q    = GA * gam
      M    = EI * chi
              
      wght = 0.5 * l0 * alpha
    
      fint  += N * bu * wght
      fint  += ( N * dot( bw , a_bar ) * bw + Q * bw ) * wght
      fint  += ( M * bt + Q * ht ) * wght
      
      stiff += EA * outer( bu , bu ) * wght
      
      stiff += EA * dot( bw , a_bar ) * outer( bu , bw ) * wght
      stiff += EA * dot( bw , a_bar ) * outer( bw , bu ) * wght
            
      stiff += ( EA * (dot( bw , a_bar ))**2 + GA + N ) * \
                 outer( bw , bw ) * wght
                 
      stiff += GA * outer( bw , ht ) * wght
      stiff += GA * outer( ht , bw ) * wght
                 
      stiff += ( EI * outer( bt , bt ) + GA * outer( ht , ht ) )* wght   
      
    stiff[4,4]=1.0                     
                                       
    elemdat.fint  = self.toGlobalCoordinates( fint  , elemdat.coords )
    elemdat.stiff = self.toGlobalCoordinates( stiff , elemdat.coords )    
    
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


  def getInternalForce ( self, elemdat ):

    EA = elemdat.props.E * elemdat.props.A
    EI = elemdat.props.E * elemdat.props.I
    l0 = norm( elemdat.coords[2]-elemdat.coords[0] )

    l0 = norm( elemdat.coords[2]-elemdat.coords[0] )
    
    intpoints = zeros(2)
    weights   = zeros(2)
    
    weights[0] = 1.
    weights[1] = 1.
       
    intpoints[0] = -0.577
    intpoints[1] =  0.577
    
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
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def getHu( self , xi ):
    
    Hu = zeros(9)
    
    Hu[0] =  0.5*(1.0-xi)
    Hu[3] = (1.0-xi*xi)
    Hu[6] =  0.5*(1.0+xi)
    
    return Hu

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def getHw( self , xi ):
    
    Hw = zeros(9)
    
    Hw[1] =  0.5*(1.0-xi)
    Hw[4] = (1.0-xi*xi)
    Hw[7] =  0.5*(1.0+xi)
    
    return Hw
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def getHt( self , xi ):
    
    Ht = zeros(9)
    
    Ht[2] =  0.5*(1.0-xi)
    Ht[5] = (1.0-xi*xi)
    Ht[8] =  0.5*(1.0+xi)
    
    return Ht
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getBu( self , l0 , xi ):

    Bu = zeros( 9 )

    Bu[0] = -1.0/l0
    Bu[3] = -4.0*xi/l0
    Bu[6] =  1.0/l0
  
    return Bu

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getBw( self , l0 , xi ):

    Bw = zeros( 9 )

    Bw[1] = -1.0/l0
    #Bw[4] = -4.0*xi/l0
    Bw[7] =  1.0/l0
    
    Bw[2] =  0.5*xi
    Bw[8] = -0.5*xi
      
    return Bw
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getBt( self , l0 , xi ):

    Bt = zeros( 9 )

    Bt[2] = -1.0/l0
    Bt[5] = -4.0*xi/l0
    Bt[8] =  1.0/l0
  
    return Bt
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def toElementCoordinates( self , a , coords ):

    a_bar = empty( a.shape )

    R     = eye( 9)
    crd   = zeros( shape=(2,2) )
    
    #crd[:,0] = coords[:,0]
    #crd[:,1] = coords[:,2]
    
    crd[0,:] = coords[0,:]
    crd[1,:] = coords[2,:]

    R[:2,:2]   = getRotationMatrix( crd )
    R[3:5,3:5] = R[:2,:2]
    R[6:8,6:8] = R[:2,:2]

    a_bar      = dot( R , a )

    return a_bar

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def toGlobalCoordinates( self , a_bar , coords ):

    a = empty( a_bar.shape )
    R     = eye( 9)
    crd   = zeros( shape=(2,2) )
    
    crd[0,:] = coords[0,:]
    crd[1,:] = coords[2,:]
    
    R[:2,:2]   = getRotationMatrix( crd )
    R[3:5,3:5] = R[:2,:2]
    R[6:8,6:8] = R[:2,:2]
    
    if len(a_bar.shape) == 1:
      a        = dot( R.transpose() , a_bar )
    elif len(a_bar.shape) == 2:
      a        = dot( R.transpose(), dot ( a_bar , R ) )
    return a

