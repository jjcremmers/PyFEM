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

from numpy import zeros, dot, eye
from scipy.linalg import norm
from math import atan2, sin, cos, tan

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class BeamNL ( Element ):

  #dofs per element
  dofTypes = [ 'u' , 'v' , 'rz' ]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.EA = self.E * self.A
    self.EI = self.E * self.I
    self.GA = self.G * self.A
    
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

    length , T = self.getT( elemdat )
        
    state = self.glob2loc( elemdat.state , T )
    
    f    = self.getF   ( state , length )
    fvar = self.getFvar( f , length )
    svar = self.getSvar( f , length )
    
    ae , kt1 = self.getTransformation( state , fvar , length ) 
    
    fint = dot( ae.transpose() , fvar )
    d    = dot( svar , ae )
    
    stiff = dot ( ae.transpose() , d )
    
    stiff = stiff + kt1
                            
    elemdat.fint  = self.loc2glob( fint  , T )
    elemdat.stiff = self.loc2glob( stiff , T )    
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):

    length , T = self.getT( elemdat )
        
    state = self.glob2loc( elemdat.state , T )
    
    f    = self.getF   ( state , length )
    fvar = self.getFvar( f , length )
    
    ae , kt1 = self.getTransformation( state , fvar , length ) 
    
    fint = dot( ae.transpose() , fvar )
                              
    elemdat.fint  = self.loc2glob( fint  , T )
         
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
    
  def getT ( self, elemdat ):

    T = zeros(shape=(2,2))
    
    lvec = elemdat.coords[1,:] - elemdat.coords[0,:]
  
    length = norm(lvec)

    lvec   = lvec*1.0/length
              
    T[0,0] =  lvec[0]
    T[0,1] =  lvec[1]
    T[1,0] = -lvec[1]
    T[1,1] =  lvec[0]
    
    return length, T
  
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getF( self , state , length ):

    
    f = zeros( 4 )

    f[0] = atan2 ( state[4] - state[1] , length + state[3] - state[0] )
    f[1] = 0.5 * ( state[2] - state[5] )
    f[2] = 0.5 * ( state[2] + state[5] ) - f[0]
    f[3] = ( state[3] - state[0] ) / length
  
    return f

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def glob2loc( self , a , t ):

    b = zeros(6)

    b[2] = a[2];
    b[5] = a[5];

    b[0] = t[0,0]*a[0] + t[0,1] * a[1];
    b[1] = a[0]*t[1,0] + t[1,1] * a[1];

    b[3] = a[3]*t[0,0] + t[0,1] * a[4];
    b[4] = a[3]*t[1,0] + t[1,1] * a[4]; 
    
    return b

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def loc2glob( self , a , T ):

    tt = eye(6)

    tt[0:2,0:2] = T
    tt[3:5,3:5] = T
    
    if a.ndim == 1:
      return dot( tt.transpose() , a )
    else:
      return dot( tt.transpose() , dot( a , tt ) )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getFvar( self , f , length ):

    fvar = zeros( 4 )

    t1 = self.EA*length
    t2 = cos(f[0])
    t3 = f[1]*f[1]
    t4 = f[2]*f[2]
    t5 = 1.0-t3/6.0-t4/10.0
    t6 = t2*t5
    t7 = 1.0+f[3]-t6
    t8 = sin(f[0])
    t12 = self.GA*length
    t13 = 1.0+f[3]
    t14 = tan(f[0])
    t17 = t13*t14-t8*t5
    t18 = t14*t14
    t25 = t7*t2
    t28 = t17*t8
    t32 = self.EI/length
    
    fvar[0] = t1*t7*t8*t5+t12*t17*(t13*(1.0+t18)-t6)
    fvar[1] = t1*t25*f[1]/3.0+t12*t28*f[1]/3.0+4.0*t32*f[1]
    fvar[2] = t1*t25*f[2]/5.0+t12*t28*f[2]/5.0+12.0*t32*f[2]
    fvar[3] = t1*t7+t12*t17*t14

    return fvar

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getSvar( self , f , length ):

    svar = zeros(shape=(4,4))
    
    t1 = self.EA*length
    t2 = sin(f[0])
    t3 = t2*t2
    t4 = f[1]*f[1]
    t5 = f[2]*f[2]
    t6 = 1.0-t4/6.0-t5/10.0
    t7 = t6*t6
    t10 = cos(f[0])
    t11 = t10*t6
    t12 = 1.0+f[3]-t11
    t13 = t12*t10
    t16 = self.GA*length
    t17 = 1.0+f[3]
    t18 = tan(f[0])
    t19 = t18*t18
    t20 = 1.0+t19
    t22 = t17*t20-t11
    t23 = t22*t22
    t25 = t17*t18
    t26 = t2*t6
    t27 = t25-t26
    t33 = t1*t2
    t36 = t12*t2
    t39 = t22*t2
    t42 = t27*t10
    t45 = t33*t11*f[1]/3.0-t1*t36*f[1]/3.0+t16*t39*f[1]/3.0+t16*t42*f[1]/3.0
    t54 = t33*t11*f[2]/5.0-t1*t36*f[2]/5.0+t16*t39*f[2]/5.0+t16*t42*f[2]/5.0
    t60 = t1*t26+t16*t22*t18+t16*t27*t20
    t61 = t10*t10
    t64 = t1*t13
    t68 = t16*t27*t2
    t70 = self.EI/length
    t78 = t1*t61*f[1]*f[2]/15.0+t16*t3*f[1]*f[2]/15.0
    t84 = t1*t10*f[1]/3.0+t16*t2*f[1]*t18/3.0
    t95 = t1*t10*f[2]/5.0+t16*t2*f[2]*t18/5.0
    
    svar[0,0] = t1*t3*t7+t1*t13*t6+t16*t23+t16*t27*(2.0*t25*t20+t26)
    svar[0,1] = t45
    svar[0,2] = t54
    svar[0,3] = t60
    svar[1,0] = t45
    svar[1,1] = t1*t61*t4/9.0+t64/3.0+t16*t3*t4/9.0+t68/3.0+4.0*t70
    svar[1,2] = t78
    svar[1,3] = t84
    svar[2,0] = t54
    svar[2,1] = t78
    svar[2,2] = t1*t61*t5/25.0+t64/5.0+t16*t3*t5/25.0+t68/5.0+12.0*t70
    svar[2,3] = t95
    svar[3,0] = t60
    svar[3,1] = t84
    svar[3,2] = t95
    svar[3,3] = t1+t16*t19
        
    return svar
    
  def getTransformation( self , u , fvar , length ):
  
    u41 = u[4]-u[1]
    Lu  = length+u[3]-u[0]

    u412 = u41*u41
    Lu2  = Lu*Lu

    u413 = u412*u41
    Lu3  = Lu2*Lu

    Lu4  = Lu2*Lu2
    Lu5  = Lu2*Lu3

    u4L  = (1.0+u412/Lu2)
    u4L2 = u4L*u4L
  
    b00  =  2.0*u41/Lu3/u4L-2.0*u413/Lu5/u4L2
    b01  = -1.0/(Lu2)/u4L+2.0/Lu4/u4L2*u412
    b11  = -2.0/Lu3/u4L2*u41

    a = zeros(shape=(4,6))
    d = zeros(shape=(6,6))

    a[0,0] =  u41/Lu2/u4L
    a[0,1] = -1.0/Lu/u4L
  
    a[0,3] = -a[0,0]
    a[0,4] = -a[0,1]

    a[1,2] =  0.5;

    a[1,5] = -0.5;
  
    a[2,0] = -a[0,0]
    a[2,1] = -a[0,1]
    a[2,2] =  0.5
  
    a[2,3] =  a[0,0]
    a[2,4] =  a[0,1]
    a[2,5] =  0.5
  
    a[3,0] = -1.0/length
    a[3,3] = -a[3,0]

  
    d[0,0] = ( fvar[0] - fvar[2] ) * b00
    d[0,1] = ( fvar[0] - fvar[2] ) * b01
    d[1,1] = ( fvar[0] - fvar[2] ) * b11

    d[0,3] = -d[0,0]
    d[0,4] = -d[0,1]
  
    d[1,0] =  d[0,1]
  
    d[1,3] = -d[0,1]
    d[1,4] = -d[1,1]
   
    d[3,0] = -d[0,0]
    d[3,1] = -d[0,1]
  
    d[3,3] =  d[0,0]
    d[3,4] =  d[0,1]
  
    d[4,0] = -d[0,1]
    d[4,1] = -d[1,1]
  
    d[4,3] =  d[0,1]
    d[4,4] =  d[1,1]
    
    return a,d

