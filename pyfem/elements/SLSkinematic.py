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

from pyfem.elements.SLSutils    import iso2loc,sigma2omega
from numpy import zeros, dot


#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class SLSkinematic:

  def __init__(self , param ):
     
    self.param    = param
 
    if param.ansFlag:
      self.epsa   = zeros( param.totDOF )
      self.epsb   = zeros( param.totDOF )
      self.epsc   = zeros( param.totDOF )
      self.epsd   = zeros( param.totDOF )
      self.epsans = zeros( shape = ( 2 , param.totDOF ) )

      self.ea     = zeros( 3 )
      self.eb     = zeros( 3 )
      self.ec     = zeros( 3 )
      self.ed     = zeros( 3 )
      self.da     = zeros( 3 )
      self.db     = zeros( 3 )
      self.dc     = zeros( 3 )
      self.dd     = zeros( 3 )

    self.du   = zeros( param.condDOF )
  
    self.dmapa13 = zeros( ( param.totDOF , param.totDOF ) )
    self.dmapb23 = zeros( ( param.totDOF , param.totDOF ) )
    self.dmapc13 = zeros( ( param.totDOF , param.totDOF ) )
    self.dmapd23 = zeros( ( param.totDOF , param.totDOF ) )
  
    class Orig:
      pass
 
    class Curr:
      pass

    class Prev:
      pass

    class Incr:
      pass

    self.orig = Orig()
    self.curr = Curr()
    self.prev = Prev()
    self.incr = Incr()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getDefVecs( self , sdat , elemdat ):
  
    '''
    
    '''
    
    psi   = sdat.h
    dpsi0 = sdat.dhdx[:,0]
    dpsi1 = sdat.dhdx[:,1]
    phi   = sdat.psi

    self.orig.x    = elemdat.coords.T
    self.orig.xBot = self.orig.x[:,:self.param.midNodes]
    self.orig.xTop = self.orig.x[:,self.param.midNodes:]

    self.curr.x    = zeros( shape=( 3 , self.param.extNodes ) )
    self.prev.x    = zeros( shape=( 3 , self.param.extNodes ) )
    self.curr.u    = zeros( shape=( 3 , self.param.extNodes ) )

    for iDim in range(3):
      for iNod in range(self.param.extNodes):
        self.curr.u[iDim,iNod] = elemdat.state[ iNod*3 + iDim ]
        
        self.curr.x[iDim,iNod] = self.orig.x[iDim,iNod] + \
                                   elemdat.state[ iNod*3 + iDim ]
        self.prev.x[iDim,iNod] = self.orig.x[iDim,iNod] + \
                                   elemdat.state0[ iNod*3 + iDim ]

    self.curr.xBot = self.curr.x[:,:self.param.midNodes]
    self.curr.xTop = self.curr.x[:,self.param.midNodes:]
    
    self.curr.uBot = self.curr.u[:,:self.param.midNodes]
    self.curr.uTop = self.curr.u[:,self.param.midNodes:]
    
    self.prev.xBot = self.prev.x[:,:self.param.midNodes]
    self.prev.xTop = self.prev.x[:,self.param.midNodes:]
   
    self.incr.xBot = self.curr.xBot - self.prev.xBot 
    self.incr.xTop = self.curr.xTop - self.prev.xTop
   
    self.curr.e1  = dot( ( self.curr.xBot + self.curr.xTop ), dpsi0 );
    self.curr.e2  = dot(   self.curr.xBot + self.curr.xTop  , dpsi1 );
    self.curr.d   = dot(  -self.curr.xBot + self.curr.xTop  , psi );
    
    self.curr.u0d1= dot( ( self.curr.uBot + self.curr.uTop ), dpsi0 );
    self.curr.u0d2= dot(   self.curr.uBot + self.curr.uTop  , dpsi1 );
    self.curr.u1  = dot(  -self.curr.uBot + self.curr.uTop  , psi );  
    self.curr.u1d1= dot(  -self.curr.uBot + self.curr.uTop  , dpsi0 );
    self.curr.u1d2= dot(  -self.curr.uBot + self.curr.uTop  , dpsi1 );
    
    self.curr.dd1 = dot(  -self.curr.xBot + self.curr.xTop  , dpsi0 );
    self.curr.dd2 = dot(  -self.curr.xBot + self.curr.xTop  , dpsi1 );

    self.curr.w   = dot(   elemdat.w , phi )

    self.prev.e1  = dot( ( self.prev.xBot + self.prev.xTop ), dpsi0 );
    self.prev.e2  = dot(   self.prev.xBot + self.prev.xTop  , dpsi1 );
    self.prev.d   = dot(  -self.prev.xBot + self.prev.xTop  , psi );
    
    self.prev.dd1 = dot(  -self.prev.xBot + self.prev.xTop , dpsi0 );
    self.prev.dd2 = dot(  -self.prev.xBot + self.prev.xTop , dpsi1 );

    self.orig.e1  = dot( ( self.orig.xBot + self.orig.xTop ), dpsi0 );
    self.orig.e2  = dot(   self.orig.xBot + self.orig.xTop  , dpsi1 );
    self.orig.d   = dot(  -self.orig.xBot + self.orig.xTop  , psi );
    
    self.orig.dd1 = dot(  -self.orig.xBot + self.orig.xTop , dpsi0 );
    self.orig.dd2 = dot(  -self.orig.xBot + self.orig.xTop , dpsi1 );
    
    self.incr.u0d1 = dot( self.incr.xBot + self.incr.xTop,  dpsi0 );
    self.incr.u0d2 = dot( self.incr.xBot + self.incr.xTop,  dpsi1 );

    self.incr.u1d1 = dot( -self.incr.xBot + self.incr.xTop,  dpsi0 );
    self.incr.u1d2 = dot( -self.incr.xBot + self.incr.xTop,  dpsi1 );

    self.incr.u1   = dot( -self.incr.xBot + self.incr.xTop,  psi );
    
    self.incr.w    = dot( elemdat.dw , phi )

    if self.param.ansFlag:
      self.getAns( sdat )
   
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def ansDmap( self , d , n1 , n2 , n3 , n4 ):
              
    for i in range(3):
      ns = (n1-1) * 3
      d[ns+i,ns+i]   = -0.125

      ns = (n2-1) * 3
      d[ns+i,ns+i]   =  0.125
    
      ns = (n3-1) * 3
      d[ns+i,ns+i]   =  0.125
    
      ns = (n4-1) * 3
      d[ns+i,ns+i]   = -0.125
    
      nis = (n1-1) * 3
      njs = (n4-1) * 3
      d[nis+i,njs+i] =  0.125
      d[njs+i,nis+i] =  0.125
    
      nis = (n2-1) * 3
      njs = (n3-1) * 3
      d[nis+i,njs+i] = -0.125
      d[njs+i,nis+i] = -0.125

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBmat( self , sdat , zeta , lamb ):
       
    bmat = zeros( ( 6 , self.param.totDOF ) )
    gbar = sdat.gbar

    psi   = sdat.h 
    dpsi0 = sdat.dhdx[:,0]
    dpsi1 = sdat.dhdx[:,1]
    phi   = sdat.psi

    for  iNod in range(self.param.midNodes):
      for iDim in range(3):
        k1=   iNod                  * 3 + iDim
        k2= ( iNod + self.param.midNodes ) * 3 + iDim
      
        bmat[0,k1] +=  self.curr.e1[iDim] * dpsi0[iNod]
        bmat[0,k2] +=  self.curr.e1[iDim] * dpsi0[iNod]
      
        bmat[1,k1] +=  self.curr.e2[iDim] * dpsi1[iNod]
        bmat[1,k2] +=  self.curr.e2[iDim] * dpsi1[iNod]
      
        bmat[2,k1] += -self.curr.d[iDim]  * psi[iNod]
        bmat[2,k2] +=  self.curr.d[iDim]  * psi[iNod]

        bmat[3,k1] +=  self.curr.e2[iDim] * dpsi0[iNod] + \
                       self.curr.e1[iDim] * dpsi1[iNod]
        bmat[3,k2] +=  self.curr.e2[iDim] * dpsi0[iNod] + \
                       self.curr.e1[iDim] * dpsi1[iNod]

        if not self.param.ansFlag: 
          bmat[4,k1] += -self.curr.e2[iDim] * psi[iNod] + \
                         self.curr.d[iDim] * dpsi1[iNod]
          bmat[4,k2] +=  self.curr.e2[iDim] * psi[iNod] + \
                         self.curr.d[iDim] * dpsi1[iNod]
      
          bmat[5,k1] += -self.curr.e1[iDim] * psi[iNod] + \
                         self.curr.d[iDim] * dpsi0[iNod]
          bmat[5,k2] +=  self.curr.e1[iDim] * psi[iNod] + \
                         self.curr.d[iDim] * dpsi0[iNod]
      
        bmat[0,k1] += zeta * \
          ( ( self.curr.dd1[iDim] - \
              2.0 * self.curr.e1[iDim] * gbar[0,0] - \
              self.curr.e2[iDim] * gbar[1,0] ) * dpsi0[iNod] - \
              self.curr.e1[iDim] * dpsi0[iNod] - \
              self.curr.e1[iDim] * gbar[1,0] * dpsi1[iNod] )
      
        bmat[0,k2] += zeta * \
          ( ( self.curr.dd1[iDim] - 2.0 * self.curr.e1[iDim] * gbar[0,0] - \
              self.curr.e2[iDim] * gbar[1,0] ) * dpsi0[iNod] + \
              self.curr.e1[iDim] * dpsi0[iNod] - \
              self.curr.e1[iDim] * gbar[1,0] * dpsi1[iNod] )
      
        bmat[1,k1] += zeta * \
          ( -self.curr.e2[iDim] * gbar[0,1] * dpsi0[iNod] + \
          ( self.curr.dd2[iDim] - 2.0 * self.curr.e2[iDim] * gbar[1,1] -
            self.curr.e1[iDim] * gbar[0,1] ) * dpsi1[iNod] -
            self.curr.e2[iDim] * dpsi1[iNod] );
      
        bmat[1,k2] += zeta * \
          ( -self.curr.e2[iDim] * gbar[0,1] * dpsi0[iNod] + \
          ( self.curr.dd2[iDim] - 2.0 * self.curr.e2[iDim] * gbar[1,1] -
            self.curr.e1[iDim] * gbar[0,1] ) * dpsi1[iNod] +
            self.curr.e2[iDim] * dpsi1[iNod] );

        bmat[2,k1] +=  zeta * 4.0 * self.curr.w * self.curr.d[iDim] * psi[iNod];

        bmat[2,k2] += -zeta * 4.0 * self.curr.w * self.curr.d[iDim] * psi[iNod];

        bmat[3,k1] +=  zeta * \
          ( ( self.curr.dd2[iDim] - self.curr.e2[iDim] * gbar[0,0] - \
              2.0 * self.curr.e1[iDim] * gbar[0,1] - \
              self.curr.e2[iDim] * gbar[1,1]) * dpsi0[iNod] - \
              self.curr.e2[iDim] * dpsi0[iNod] + \
            ( self.curr.dd1[iDim] - self.curr.e1[iDim] * gbar[0,0] - \
              2.0 * self.curr.e2[iDim] * gbar[1,0] - \
              self.curr.e1[iDim] * gbar[1,1]) * dpsi1[iNod] - \
              self.curr.e1[iDim] * dpsi1[iNod] ) 

        bmat[3,k2] += zeta * \
          ( ( self.curr.dd2[iDim] - self.curr.e2[iDim] * gbar[0,0] - \
             2.0 * self.curr.e1[iDim] * gbar[0,1] - \
              self.curr.e2[iDim] * gbar[1,1]) * dpsi0[iNod] + \
              self.curr.e2[iDim] * dpsi0[iNod] + \
            ( self.curr.dd1[iDim] - self.curr.e1[iDim] * gbar[0,0] - \
              2.0 * self.curr.e2[iDim] * gbar[1,0] - \
              self.curr.e1[iDim] * gbar[1,1]) * dpsi1[iNod] +
              self.curr.e1[iDim] * dpsi1[iNod] );

        bmat[4,k1] += zeta * (-self.curr.dd2[iDim] * psi[iNod,] - \
                               self.curr.d[iDim] * dpsi1[iNod] )

        bmat[4,k2] += zeta * ( self.curr.dd2[iDim] * psi[iNod] + \
                               self.curr.d[iDim] * dpsi1[iNod] )

        bmat[5,k1] += zeta * (-self.curr.dd1[iDim] * psi[iNod] - \
                               self.curr.d[iDim] * dpsi0[iNod] )

        bmat[5,k2] += zeta * ( self.curr.dd1[iDim] * psi[iNod] + \
                               self.curr.d[iDim] * dpsi0[iNod] )
   
    dnorm = dot( self.curr.d , self.curr.d )

    for k1 in range(self.param.intNodes):
      bmat[ 2, self.param.condDOF+k1 ] += -2.0 * zeta * dnorm * phi[k1];
  
    #Assumed natural strains                  

    if self.param.ansFlag:
      bmat[4,:] = self.epsans[0,:]
      bmat[5,:] = self.epsans[1,:]

    return iso2loc( bmat , lamb )
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getStrains( self , kin , sdat , zeta , lamb ):

    kin.strain  = iso2loc( self.getEps()  + zeta * self.getRho ( sdat.gbar ) , lamb )
    kin.dstrain = iso2loc( self.getDEps() + zeta * self.getDRho( sdat.gbar ) , lamb )
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
  def getEps( self ):

    eps = zeros(6)

    eps[0] = 0.5 * ( dot( self.orig.e1   , self.curr.u0d1 ) + \
                     dot( self.orig.e1   , self.curr.u0d1 ) + \
                     dot( self.curr.u0d1 , self.curr.u0d1 ) )

    eps[1] = 0.5 * ( dot( self.orig.e2   , self.curr.u0d2 ) + \
                     dot( self.orig.e2   , self.curr.u0d2 ) + \
                     dot( self.curr.u0d2 , self.curr.u0d2 ) )

    eps[2] =         dot( self.orig.d    , self.curr.u1   ) + \
               0.5 * dot( self.curr.u1   , self.curr.u1   )
    

    eps[3] =         dot( self.orig.e2   , self.curr.u0d1 ) + \
                     dot( self.orig.e1   , self.curr.u0d2 ) + \
                     dot( self.curr.u0d1 , self.curr.u0d2 )

    eps[4] =         dot( self.orig.e2   , self.curr.u1   ) + \
                     dot( self.orig.d    , self.curr.u0d2 ) + \
                     dot( self.curr.u0d2 , self.curr.u1   )

    eps[5] =         dot( self.orig.e1   , self.curr.u1   ) + \
                     dot( self.orig.d    , self.curr.u0d1 ) + \
                     dot( self.curr.u0d1 , self.curr.u1   )

    return eps
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
  def getDEps( self ):

    deps = zeros(6)

    deps[0] = 0.5 * ( dot( self.prev.e1   , self.incr.u0d1 ) + \
                      dot( self.prev.e1   , self.incr.u0d1 ) + \
                      dot( self.incr.u0d1 , self.incr.u0d1 ) )

    deps[1] = 0.5 * ( dot( self.prev.e2   , self.incr.u0d2 ) + \
                      dot( self.prev.e2   , self.incr.u0d2 ) + \
                      dot( self.incr.u0d2 , self.incr.u0d2 ) )

    deps[2] =         dot( self.prev.d    , self.incr.u1   ) + \
                0.5 * dot( self.incr.u1   , self.incr.u1   )

    deps[3] =         dot( self.prev.e2   , self.incr.u0d1 ) + \
                      dot( self.prev.e1   , self.incr.u0d2 ) + \
                      dot( self.incr.u0d1 , self.incr.u0d2 )

    deps[4] =         dot( self.prev.e2   , self.incr.u1   ) + \
                      dot( self.prev.d    , self.incr.u0d2 ) + \
                      dot( self.incr.u0d2 , self.incr.u1   )

    deps[5] =         dot( self.prev.e1   , self.incr.u1   ) + \
                      dot( self.prev.d    , self.incr.u0d1 ) + \
                      dot( self.incr.u0d1 , self.incr.u1   )

    return deps


#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getRho( self , gbar ):
 
    rho = zeros(6)
    
    rho[0] = 0.5 * ( dot( self.orig.dd1  , self.curr.u0d1 ) + \
                     dot( self.curr.u1d1 , self.orig.e1   ) + \
                     dot( self.curr.u1d1 , self.curr.u0d1 ) + \
                     dot( self.orig.dd1  , self.curr.u0d1 ) + \
                     dot( self.curr.u1d1 , self.orig.e1   ) + \
                     dot( self.curr.u1d1 , self.curr.u0d1 ) )
 
    rho[1] = 0.5 * ( dot( self.orig.dd2   , self.curr.u0d2 ) + \
                     dot( self.curr.u1d2  , self.orig.e2   ) + \
                     dot( self.orig.dd2   , self.curr.u0d2 ) + \
                     dot( self.curr.u1d2  , self.orig.e2   ) + \
                     dot( self.curr.u1d2  , self.curr.u0d2 ) )

    rho[2] = -2.0 * dot( self.curr.d , self.curr.d ) * self.curr.w
 
    rho[3] =        dot( self.orig.dd1  , self.curr.u0d2 ) + \
                    dot( self.curr.u1d1 , self.orig.e2   ) + \
                    dot( self.curr.u1d1 , self.curr.u0d2 ) + \
                    dot( self.orig.dd2  , self.curr.u0d1 ) + \
                    dot( self.curr.u1d2 , self.orig.e1   ) + \
                    dot( self.curr.u1d2 , self.curr.u0d1 )


    rho[4] =        dot( self.orig.dd2  , self.curr.u1 ) + \
                    dot( self.curr.u1d2 , self.orig.d  ) + \
                    dot( self.curr.u1d2 , self.curr.u1 )

    rho[5] =        dot( self.orig.dd1  , self.curr.u1 ) + \
                    dot( self.curr.u1d1 , self.orig.d  ) + \
                    dot( self.curr.u1d1 , self.curr.u1 )
    
    return rho
	
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getDRho( self , gbar ):
 
    drho = zeros(6)
    
    drho[0] = 0.5 * ( dot ( self.prev.e1   , self.incr.u1d1 ) + \
                      dot ( self.prev.dd1  , self.incr.u0d1 ) + \
                      dot ( self.incr.u1d1 , self.incr.u0d1 ) + \
                      dot ( self.prev.e1   , self.incr.u1d1 ) + \
                      dot ( self.prev.dd1  , self.incr.u0d1 ) + \
                      dot ( self.incr.u1d1 , self.incr.u0d1 ) )
    
    drho[1] = 0.5 * ( dot ( self.prev.e2   , self.incr.u1d2 ) + \
                      dot ( self.prev.dd2  , self.incr.u0d2 ) + \
                      dot ( self.incr.u1d2 , self.incr.u0d2 ) + \
                      dot ( self.prev.e2   , self.incr.u1d2 ) + \
                      dot ( self.prev.dd2  , self.incr.u0d2 ) + \
                      dot ( self.incr.u1d2 , self.incr.u0d2 ) )

    drho[2] = -4.0 * self.curr.w * dot ( self.prev.d  , self.incr.u1 ) + \
              -2.0 * self.curr.w * dot ( self.incr.u1 , self.incr.u1 ) + \
              -2.0 * self.incr.w * dot ( self.prev.d  , self.prev.d  ) + \
              -4.0 * self.incr.w * dot ( self.incr.u1 , self.prev.d  ) + \
              -2.0 * self.incr.w * dot ( self.incr.u1 , self.incr.u1 )
  
    drho[3] = dot ( self.prev.e2   , self.incr.u1d1 ) + \
              dot ( self.prev.dd1  , self.incr.u0d2 ) + \
              dot ( self.incr.u1d1 , self.incr.u0d2 ) + \
              dot ( self.prev.e1   , self.incr.u1d2 ) + \
              dot ( self.prev.dd2  , self.incr.u0d1 ) + \
              dot ( self.incr.u1d2 , self.incr.u0d1 )
  
    drho[4] = dot ( self.prev.dd2  , self.incr.u1   ) + \
              dot ( self.prev.d    , self.incr.u1d2 ) + \
              dot ( self.incr.u1d2 , self.incr.u1   )

    drho[5] = dot ( self.prev.dd1  , self.incr.u1   ) + \
              dot ( self.prev.d    , self.incr.u1d1 ) + \
              dot ( self.incr.u1d1 , self.incr.u1   )
  
    drho[0] += - gbar[0,0] * ( dot( self.prev.e1   , self.incr.u0d1 ) + \
                               dot( self.prev.e1   , self.incr.u0d1 ) + \
                               dot( self.incr.u0d1 , self.incr.u0d1 ) ) - \
                 gbar[1,0] * ( dot( self.prev.e2   , self.incr.u0d1 ) + \
                               dot( self.prev.e1   , self.incr.u0d2 ) + \
                               dot( self.incr.u0d2 , self.incr.u0d1 ) ) - \
                 gbar[0,0] * ( dot( self.prev.e1   , self.incr.u0d1 ) + \
                               dot( self.prev.e1   , self.incr.u0d1 ) + \
                               dot( self.incr.u0d1 , self.incr.u0d1 ) ) - \
                 gbar[1,0] * ( dot( self.prev.e1   , self.incr.u0d2 ) + \
                               dot( self.prev.e2   , self.incr.u0d1 ) + \
                               dot( self.incr.u0d1 , self.incr.u0d2 ) )

    drho[1] += - 0.5*gbar[0,1] * ( dot( self.prev.e1   , self.incr.u0d2 ) + \
                                   dot( self.prev.e2   , self.incr.u0d1 ) + \
                                   dot( self.incr.u0d1 , self.incr.u0d2 ) ) - \
                0.5*gbar[1,1] * ( dot( self.prev.e2   , self.incr.u0d2 ) + \
                              dot( self.prev.e2   , self.incr.u0d2 ) + \
                              dot( self.incr.u0d2 , self.incr.u0d2 ) ) - \
                0.5*gbar[0,1] * ( dot( self.prev.e2   , self.incr.u0d1 ) + \
                              dot( self.prev.e1   , self.incr.u0d2 ) + \
                              dot( self.incr.u0d1 , self.incr.u0d2 ) ) - \
                0.5*gbar[1,1] * ( dot( self.prev.e2   , self.incr.u0d2 ) + \
                              dot( self.prev.e2   , self.incr.u0d2 ) + \
                              dot( self.incr.u0d2 , self.incr.u0d2 ) )

    drho[3] += - gbar[0,0] * ( dot( self.prev.e1   , self.incr.u0d2 ) + \
                              dot( self.prev.e2   , self.incr.u0d1 ) + \
                              dot( self.incr.u0d1 , self.incr.u0d2 ) ) - \
                gbar[1,0] * ( dot( self.prev.e2   , self.incr.u0d2 ) + \
                              dot( self.prev.e2   , self.incr.u0d2 ) + \
                              dot( self.incr.u0d2 , self.incr.u0d2 ) ) - \
                gbar[0,1] * ( dot( self.prev.e1   , self.incr.u0d1 ) + \
                              dot( self.prev.e1   , self.incr.u0d1 ) + \
                              dot( self.incr.u0d1 , self.incr.u0d1 ) ) - \
                gbar[1,1] * ( dot( self.prev.e1   , self.incr.u0d2 ) + \
                              dot( self.prev.e2   , self.incr.u0d1 ) + \
                              dot( self.incr.u0d1 , self.incr.u0d2 ) )

    return drho

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getAns( self , sdat ):
        
    fa = 0.5 * ( 1.0 - sdat.xi[1] )
    fb = 0.5 * ( 1.0 + sdat.xi[0] )
    fc = 0.5 * ( 1.0 + sdat.xi[1] )
    fd = 0.5 * ( 1.0 - sdat.xi[0] )
  
    for i in range(3):
      self.ea[i] = -self.curr.x[i,4]-self.curr.x[i,0]+self.curr.x[i,5]+self.curr.x[i,1]
      self.eb[i] = -self.curr.x[i,5]-self.curr.x[i,1]+self.curr.x[i,6]+self.curr.x[i,2]
      self.ec[i] = -self.curr.x[i,7]-self.curr.x[i,3]+self.curr.x[i,6]+self.curr.x[i,2]
      self.ed[i] = -self.curr.x[i,4]-self.curr.x[i,0]+self.curr.x[i,7]+self.curr.x[i,3]
    
      self.da[i] =  self.curr.x[i,4]-self.curr.x[i,0]+self.curr.x[i,5]-self.curr.x[i,1]
      self.db[i] =  self.curr.x[i,5]-self.curr.x[i,1]+self.curr.x[i,6]-self.curr.x[i,2]
      self.dc[i] =  self.curr.x[i,7]-self.curr.x[i,3]+self.curr.x[i,6]-self.curr.x[i,2]
      self.dd[i] =  self.curr.x[i,4]-self.curr.x[i,0]+self.curr.x[i,7]-self.curr.x[i,3]

    for iDim in range(3):
      self.epsa[4 * 3 + iDim] =  self.ea[iDim] - self.da[iDim]
      self.epsa[0 * 3 + iDim] = -self.ea[iDim] - self.da[iDim] 
      self.epsa[5 * 3 + iDim] =  self.ea[iDim] + self.da[iDim]     
      self.epsa[1 * 3 + iDim] = -self.ea[iDim] + self.da[iDim]
    
      self.epsb[5 * 3 + iDim] =  self.eb[iDim] - self.db[iDim]
      self.epsb[1 * 3 + iDim] = -self.eb[iDim] - self.db[iDim]
      self.epsb[6 * 3 + iDim] =  self.eb[iDim] + self.db[iDim]
      self.epsb[2 * 3 + iDim] = -self.eb[iDim] + self.db[iDim]
    
      self.epsc[7 * 3 + iDim] =  self.ec[iDim] - self.dc[iDim]
      self.epsc[3 * 3 + iDim] = -self.ec[iDim] - self.dc[iDim]
      self.epsc[6 * 3 + iDim] =  self.ec[iDim] + self.dc[iDim]
      self.epsc[2 * 3 + iDim] = -self.ec[iDim] + self.dc[iDim]
    
      self.epsd[4 * 3 + iDim] =  self.ed[iDim] - self.dd[iDim]
      self.epsd[0 * 3 + iDim] = -self.ed[iDim] - self.dd[iDim]
      self.epsd[7 * 3 + iDim] =  self.ed[iDim] + self.dd[iDim]
      self.epsd[3 * 3 + iDim] = -self.ed[iDim] + self.dd[iDim]
  
    self.epsans[0,:] = 0.0625 * ( fb * self.epsb + fd * self.epsd )
    self.epsans[1,:] = 0.0625 * ( fa * self.epsa + fc * self.epsc )
  
    self.ansDmap( self.dmapa13 , 5 , 1 , 6 , 2 )
    self.ansDmap( self.dmapb23 , 6 , 2 , 7 , 3 )
    self.ansDmap( self.dmapc13 , 8 , 4 , 7 , 3 )
    self.ansDmap( self.dmapd23 , 5 , 1 , 8 , 4 )
  
    self.d13 = fa * self.dmapa13 + fc * self.dmapc13
    self.d23 = fd * self.dmapd23 + fb * self.dmapb23

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def addGeomStiff( self , stiff , sdat , sigma , lamb, z ):   
    
    omega = sigma2omega( sigma , lamb )

    for iNod in range( self.param.extNodes ):
      i = 3*iNod
   
      if iNod < self.param.midNodes:
        pi1i   = -sdat.h[iNod]
        pi0d1i =  sdat.dhdx[iNod,0]
        pi0d2i =  sdat.dhdx[iNod,1]
        pi1d1i = -sdat.dhdx[iNod,0]
        pi1d2i = -sdat.dhdx[iNod,1]
      else:
        pi1i   =  sdat.h[iNod-self.param.midNodes]
        pi0d1i =  sdat.dhdx[iNod-self.param.midNodes,0]
        pi0d2i =  sdat.dhdx[iNod-self.param.midNodes,1]
        pi1d1i =  sdat.dhdx[iNod-self.param.midNodes,0]
        pi1d2i =  sdat.dhdx[iNod-self.param.midNodes,1]

      for jNod in range( self.param.extNodes ):
        j = 3*jNod

        if jNod < self.param.midNodes:
          pi1j   = -sdat.h[jNod]
          pi0d1j =  sdat.dhdx[jNod,0]
          pi0d2j =  sdat.dhdx[jNod,1]
          pi1d1j = -sdat.dhdx[jNod,0]
          pi1d2j = -sdat.dhdx[jNod,1]
        else:
          pi1j   =  sdat.h[jNod-self.param.midNodes]
          pi0d1j =  sdat.dhdx[jNod-self.param.midNodes,0]
          pi0d2j =  sdat.dhdx[jNod-self.param.midNodes,1]
          pi1d1j =  sdat.dhdx[jNod-self.param.midNodes,0]
          pi1d2j =  sdat.dhdx[jNod-self.param.midNodes,1]

        add  = omega[0] * pi0d1i * pi0d1j

        add += omega[1] * pi0d2i * pi0d2j
  
        add += omega[2] * pi1i * pi1j

        add += omega[3] * (pi0d1i * pi0d2j + pi0d1j * pi0d2i)

        if not self.param.ansFlag:
          add += omega[4] * (pi0d2i * pi1j + pi0d2j * pi1i)        
          add += omega[5] * (pi0d1i * pi1j + pi0d1j * pi1i)
      
        add +=  z * omega[0] * (pi1d1i * pi0d1j + pi1d1j * pi0d1i)
        add += -z * sdat.gbar[0,0] * omega[0] * (pi0d1i * pi0d1j + pi0d1j * pi0d1i)
        add += -z * sdat.gbar[1,0] * omega[0] * (pi0d1i * pi0d2j + pi0d1j * pi0d2i)
  
        add +=  z * omega[1] * (pi1d2i * pi0d2j + pi1d2j * pi0d2i)
        add += -z * sdat.gbar[0,1] * omega[1] * (pi0d1i * pi0d2j + pi0d1j * pi0d2i)
        add += -z * sdat.gbar[1,1] * omega[1] * (pi0d2i * pi0d2j + pi0d2j * pi0d2i)
  
        add += -z * 4.0 * self.curr.w * omega[2] * (pi1i * pi1j)
    
        add +=  z * omega[4] * (pi1d2i * pi1j + pi1d2j * pi1i)
  
        add +=  z * omega[5] * (pi1d1i * pi1j + pi1d1j * pi1i)
   
        add +=  z * omega[3] * (pi1d1i * pi0d2j + pi1d1j * pi0d2i)
        add +=  z * omega[3] * (pi1d2i * pi0d1j + pi1d2j * pi0d1i)
        add += -z * (sdat.gbar[0,0] + sdat.gbar[1,1]) * \
                     omega[3] * (pi0d1i * pi0d2j + pi0d1j * pi0d2i)
        add += -z * sdat.gbar[1,0] * omega[3] * (pi0d2i * pi0d2j + pi0d2j * pi0d2i)
        add += -z * sdat.gbar[0,1] * omega[3] * (pi0d1i * pi0d1j + pi0d1j * pi0d1i)

        stiff[i+0,j+0] += add
        stiff[i+1,j+1] += add
        stiff[i+2,j+2] += add
  
    fac = -z * 4.0 * omega[2]
   
    self.geom04( sdat , self.curr.d , fac , stiff );
 
    if self.param.ansFlag:
      stiff += omega[4] * self.d23 + omega[5] * self.d13

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def geom04( self , sdat , d , fac , stiff ):
       
    for i,psi in enumerate(sdat.psi):
      tmp = psi * fac
      ii  = i + self.param.condDOF

      for j in range(self.param.extNodes):
        if j < self.param.midNodes:
          pi1 = -sdat.h[j]    
        else:
          pi1 =  sdat.h[j-self.param.midNodes]
          
        add = tmp * pi1
        jj  = j   * self.param.nodeDOF

        stiff[ii,jj+0] += d[0] * add;
        stiff[ii,jj+1] += d[1] * add;
        stiff[ii,jj+2] += d[2] * add;
      
        stiff[jj+0,ii] += d[0] * add;
        stiff[jj+1,ii] += d[1] * add;
        stiff[jj+2,ii] += d[2] * add;

