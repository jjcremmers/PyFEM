from numpy import dot,zeros,insert,array,ones
from math import sqrt

def vonMisesStress( s ):

  smises = 0.;

  if len(s) == 3:
    return sqrt(s[0]*s[0]+s[1]*s[1]-s[0]*s[1]+3.*s[2]*s[2]);
  elif len(s) == 6:
    smises = ( s[0] - s[1] ) * ( s[0] - s[1] ) + \
             ( s[1] - s[2] ) * ( s[1] - s[2] ) + \
             ( s[2] - s[0] ) * ( s[2] - s[0] )
 
    smises += 6.0 * dot(s[3:] , s[3:] )
    return sqrt( 0.5 * smises )

def hydrostaticStress( s ):

  return 0.333333333333333*sum( s[:3] );

class Hardening:

  def __init__ ( self, props ):

    self.n = 20
    self.maxStrain = 1.0

    for name,val in props:
      setattr( self, name, val )

    self.htype = 0

    if hasattr( props , "EqPlasStrains" ):
      self.htype = 1
      insert(self.Stresses, 0, self.syield )
      insert(self.EqPlasStrains, 0, 0. )
      
    elif hasattr( props , "q" ):
      self.htype = 2
      self.Stresses = zeros( self.n+1 )
      self.EqPlasStrains = zeros( self.n + 1 )
    
      epsInc = self.maxStrain / ( 1.0*self.n )

      if not hasattr( props , "K" ):
        self.K = self.syield*pow(self.E/self.syield,self.q)

      eps0 = self.syield/self.E

      self.Stresses[0] = self.syield

      for i in range(self.n):
        self.EqPlasStrains[i+1] = (i+1)*epsInc
        self.Stresses[i+1] = self.K * pow(self.EqPlasStrains[i+1]+eps0,self.q)

    print("RR",self.Stresses,self.EqPlasStrains)
 
  def getHardening( self , eqplas ):

#    if self.n == 0:
#      return self.Stresses[0] , self.EqPlasStrains[0]
#    else:

    
    if True:

      for i in range(self.n):
        eqpl1 = self.EqPlasStrains[i+1]

        if eqplas < eqpl1:
          eqpl0  = self.EqPlasStrains[i]
          deqpl  = eqpl1-eqpl0
    
          syiel0 = self.Stresses[i]
          dsyiel = self.Stresses[i+1]-syiel0
          
          hard = dsyiel / deqpl
          return syiel0 + ( eqplas-eqpl0)*hard  , hard
     
      eqpl0  = self.EqPlasStrains[:-1]
      deqpl  = self.EqPlasStrains[:]-eqpl0
    
      syiel0 = self.Stresses[:-1]
      dsyiel = self.Stresses[:]-syiel0
     
      hard = dsyiel / deqpl
      return siyel0 + ( eqplas-eqpl0)*hard  , hard
      
#
#
#

def transform2To3( s ):
  return array([ s[0] , s[1] , 0. , 0. , 0. , s[2] ])

def transform3To2( s , t ):
  return array([ s[0] , s[1] , s[5] ]) , \
         array([(t[0,0] , t[0,1] , t[0,5]), \
                (t[1,0] , t[1,1] , t[1,5]), \
                (t[5,0] , t[5,1] , t[5,5])])
                
                
#------------------------------------------------------------------------------
#  hydroStaticStress
#------------------------------------------------------------------------------


def hydroStatic( s ):

  return sum( s[:3] )/3.0


    
#------------------------------------------------------------------------------
#  splitStrains
#------------------------------------------------------------------------------


def splitStrains( eps ):

  theta = hydroStatic( eps )
  
  deps = zeros(6)

  deps[:3] = eps[:3] - theta*ones(3)
  deps[3:] = 0.5*eps[3:]

  return theta,deps  
                   

    
