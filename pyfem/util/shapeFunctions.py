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
from math import sqrt
from numpy import array, dot, ndarray, empty, zeros , ones
from scipy.linalg import norm , det , inv
from scipy.special.orthogonal import p_roots as gauss_scheme

class shapeData:
  
  pass
   
#----------------------------------------------------------------------

class elemShapeData:

  def __init__( self ):
    
    self.sData = []

  def __iter__( self ):

    return iter(self.sData)

  def __len__( self ):

    return len(self.sData)
            
#----------------------------------------------------------------------

def getShapeLine2 ( xi ):

  #Check the dimensions of the physical space
  if type(xi) != float:
    raise NotImplementedError('1D only')

  sData       = shapeData()
  
  #Set length of lists
  sData.h     = empty( 2 )
  sData.dhdxi = empty( shape=(2,1) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = 0.5*(1.0-xi)
  sData.h[1] = 0.5*(1.0+xi)

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.5
  sData.dhdxi[1,0] =  0.5

  return sData

#----------------------------------------------------------------------

def getShapeLine3 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 1:
    raise NotImplementedError('1D only')

  sData       = shapeData()
  
  #Set length of lists
  sData.h     = empty( 3 )
  sData.dhdxi = empty( shape=(1,3) )
  sData.xi    = xi

  #Calculate shape functions
  phi[0] = 0.5*(1.0-xi)-0.5*(1.0-xi*xi)
  phi[1] = 1-xi*xi
  phi[2] = 0.5*(1.0+xi)-0.5*(1.0-xi*xi)

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.5+xi
  sData.dhdxi[0,1] = -2.0*xi
  sData.dhdxi[0,2] =  0.5+xi

  return sData

#----------------------------------------------------------------------

def getShapeTria3 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 2:
    raise NotImplementedError('2D only')

  sData       = shapeData()
  
  #Set length of lists
  sData.h     = empty( 3 )
  sData.dhdxi = empty( shape=(3,2) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = 1.0-xi[0]-xi[1]
  sData.h[1] = xi[0]
  sData.h[2] = xi[1]

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -1.0
  sData.dhdxi[1,0] =  1.0
  sData.dhdxi[2,0] =  0.0

  sData.dhdxi[0,1] = -1.0
  sData.dhdxi[1,1] =  0.0
  sData.dhdxi[2,1] =  1.0

  return sData

#-------------------------------------

def getShapeQuad4 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 2:
    raise NotImplementedError('2D only')

  sData       = shapeData()
  
  #Set length of lists
  sData.h     = empty( 4 )
  sData.dhdxi = empty( shape=(4,2) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = 0.25*(1.0-xi[0])*(1.0-xi[1])
  sData.h[1] = 0.25*(1.0+xi[0])*(1.0-xi[1])
  sData.h[2] = 0.25*(1.0+xi[0])*(1.0+xi[1])
  sData.h[3] = 0.25*(1.0-xi[0])*(1.0+xi[1])

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.25*(1.0-xi[1])
  sData.dhdxi[1,0] =  0.25*(1.0-xi[1])
  sData.dhdxi[2,0] =  0.25*(1.0+xi[1])
  sData.dhdxi[3,0] = -0.25*(1.0+xi[1])

  sData.dhdxi[0,1] = -0.25*(1.0-xi[0])
  sData.dhdxi[1,1] = -0.25*(1.0+xi[0])
  sData.dhdxi[2,1] =  0.25*(1.0+xi[0])
  sData.dhdxi[3,1] =  0.25*(1.0-xi[0])

  return sData

#-------------------------------------

def getShapeTria6 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 2:
    raise NotImplementedError('2D only')

  sData       = shapeData()
   
  #Set length of lists
  sData.h     = empty( 6 )
  sData.dhdxi = empty( shape=(6,2) )
  sData.xi    = xi

  sData.h[0] = 1.0-xi[0]-xi[1]
  sData.h[1] = xi[0]
  sData.h[2] = xi[1]

  #Calculate shape functions
  sData.h[0] = 1.0-xi[0]-xi[1]-2.0*xi[0]*(1.0-xi[0]-xi[1])-2.0*xi[1]*(1.0-xi[0]-xi[1])
  sData.h[1] = xi[0]-2.0*xi[0]*(1.0-xi[0]-xi[1])-2.0*xi[0]*xi[1]
  sData.h[2] = xi[1]-2.0*xi[0]*xi[1]-2.0*xi[1]*(1.0-xi[0]-xi[1])
  sData.h[3] = 4.0*xi[0]*(1.0-xi[0]-xi[1])
  sData.h[4] = 4.0*xi[0]*xi[1]
  sData.h[5] = 4.0*xi[1]*(1.0-xi[0]-xi[1])

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -1.0-2.0*(1.0-xi[0]-xi[1])+2.0*xi[0]+2.0*xi[1]
  sData.dhdxi[1,0] =  1.0-2.0*(1.0-xi[0]-xi[1])+2.0*xi[0]-2.0*xi[1]
  sData.dhdxi[2,0] =  0.0
  sData.dhdxi[3,0] =  4.0*(1.0-xi[0]-xi[1])-4.0*xi[0]
  sData.dhdxi[4,0] =  4.0*xi[1]
  sData.dhdxi[5,0] = -4.0*xi[1]

  sData.dhdxi[0,1] = -1.0+2.0*xi[0]-2.0*(1.0-xi[0]-xi[1])+2.0*xi[1]
  sData.dhdxi[1,1] =  0.0
  sData.dhdxi[2,1] =  1.0-2.0*xi[0]-2.0*(1.0-xi[0]-xi[1])+2.0*xi[1]
  sData.dhdxi[3,1] = -4.0*xi[0]
  sData.dhdxi[4,1] =  4.0*xi[0]
  sData.dhdxi[5,1] =  4.0*(1.0-xi[0]-xi[1])-4.0*xi[1]

  return sData

#-------------------------------------

def getShapeQuad8 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 2:
    raise NotImplementedError('2D only')

  sData       = shapeData()
  
  #Set length of lists
  sData.h     = empty( 8 )
  sData.dhdxi = empty( shape=(8,2) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = -0.25*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[0]+xi[1])
  sData.h[1] =  0.5 *(1.0-xi[0])*(1.0+xi[0])*(1.0-xi[1])
  sData.h[2] = -0.25*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[0]+xi[1])
  sData.h[3] =  0.5 *(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[1])
  sData.h[4] = -0.25*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[0]-xi[1])
  sData.h[5] =  0.5 *(1.0-xi[0])*(1.0+xi[0])*(1.0+xi[1])
  sData.h[6] = -0.25*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[0]-xi[1])
  sData.h[7] =  0.5 *(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[1])

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.25*(-1.0+xi[1])*( 2.0*xi[0]+xi[1])
  sData.dhdxi[1,0] =  xi[0]*(-1.0+xi[1])
  sData.dhdxi[2,0] =  0.25*(-1.0+xi[1])*(-2.0*xi[0]+xi[1])
  sData.dhdxi[3,0] = -0.5 *(1.0+xi[1])*(-1.0+xi[1])
  sData.dhdxi[4,0] =  0.25*( 1.0+xi[1])*( 2.0*xi[0]+xi[1])
  sData.dhdxi[5,0] = -xi[0]*(1.0+xi[1])
  sData.dhdxi[6,0] = -0.25*( 1.0+xi[1])*(-2.0*xi[0]+xi[1])
  sData.dhdxi[7,0] = 0.5*(1.0+xi[1])*(-1.0+xi[1])

  sData.dhdxi[0,1] = -0.25*(-1.0+xi[0])*( xi[0]+2.0*xi[1])
  sData.dhdxi[1,1] =  0.5 *( 1.0+xi[0])*(-1.0+xi[0])
  sData.dhdxi[2,1] =  0.25*( 1.0+xi[0])*(-xi[0]+2.0*xi[1])
  sData.dhdxi[3,1] = -xi[1]*(1.0+xi[0])
  sData.dhdxi[4,1] =  0.25*( 1.0+xi[0])*( xi[0]+2.0*xi[1])
  sData.dhdxi[5,1] = -0.5 *( 1.0+xi[0])*(-1.0+xi[0])
  sData.dhdxi[6,1] = -0.25*(-1.0+xi[0])*(-xi[0]+2.0*xi[1])
  sData.dhdxi[7,1] =  xi[1]*(-1.0+xi[0])

  return sData

#-------------------------------------

def getShapeQuad9 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 2:
    raise NotImplementedError('2D only')

  sData       = shapeData()
  
  #Set length of lists
  sData.h     = empty( 9 )
  sData.dhdxi = empty( shape=(9,2) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = -0.25*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[0]+xi[1])
  sData.h[1] =  0.5 *(1.0-xi[0])*(1.0+xi[0])*(1.0-xi[1])
  sData.h[2] = -0.25*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[0]+xi[1])
  sData.h[3] =  0.5 *(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[1])
  sData.h[4] = -0.25*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[0]-xi[1])
  sData.h[5] =  0.5 *(1.0-xi[0])*(1.0+xi[0])*(1.0+xi[1])
  sData.h[6] = -0.25*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[0]-xi[1])
  sData.h[7] =  0.5 *(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[1])
  sData.h[8] =  0.5 *(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[1])

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.25*(-1.0+xi[1])*( 2.0*xi[0]+xi[1])
  sData.dhdxi[1,0] =  xi[0]*(-1.0+xi[1])
  sData.dhdxi[2,0] =  0.25*(-1.0+xi[1])*(-2.0*xi[0]+xi[1])
  sData.dhdxi[3,0] = -0.5 *(1.0+xi[1])*(-1.0+xi[1])
  sData.dhdxi[4,0] =  0.25*( 1.0+xi[1])*( 2.0*xi[0]+xi[1])
  sData.dhdxi[5,0] = -xi[0]*(1.0+xi[1])
  sData.dhdxi[6,0] = -0.25*( 1.0+xi[1])*(-2.0*xi[0]+xi[1])
  sData.dhdxi[7,0] = 0.5*(1.0+xi[1])*(-1.0+xi[1])
  sData.dhdxi[8,0] = 0.5*(1.0+xi[1])*(-1.0+xi[1])

  sData.dhdxi[0,1] = -0.25*(-1.0+xi[0])*( xi[0]+2.0*xi[1])
  sData.dhdxi[1,1] =  0.5 *( 1.0+xi[0])*(-1.0+xi[0])
  sData.dhdxi[2,1] =  0.25*( 1.0+xi[0])*(-xi[0]+2.0*xi[1])
  sData.dhdxi[3,1] = -xi[1]*(1.0+xi[0])
  sData.dhdxi[4,1] =  0.25*( 1.0+xi[0])*( xi[0]+2.0*xi[1])
  sData.dhdxi[5,1] = -0.5 *( 1.0+xi[0])*(-1.0+xi[0])
  sData.dhdxi[6,1] = -0.25*(-1.0+xi[0])*(-xi[0]+2.0*xi[1])
  sData.dhdxi[7,1] =  xi[1]*(-1.0+xi[0])
  sData.dhdxi[8,1] =  xi[1]*(-1.0+xi[0])
  return sData

#----------------------------------------------------------------------

def getShapeTetra4 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 3:
    raise NotImplementedError('3D only')

  sData       = shapeData()
 
  #Set length of lists
  sData.h     = empty( 4 )
  sData.dhdxi = empty( shape=(4,3) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = 0.25*(1.0-xi[0])*(1.0-xi[1])
  sData.h[1] = 0.25*(1.0+xi[0])*(1.0-xi[1])
  sData.h[2] = 0.25*(1.0+xi[0])*(1.0+xi[1])
  sData.h[3] = 0.25*(1.0-xi[0])*(1.0+xi[1])

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.25*(1.0-xi[1])
  sData.dhdxi[1,0] =  0.25*(1.0-xi[1])
  sData.dhdxi[2,0] =  0.25*(1.0+xi[1])
  sData.dhdxi[3,0] = -0.25*(1.0+xi[1])

  sData.dhdxi[0,1] = -0.25*(1.0-xi[0])
  sData.dhdxi[1,1] = -0.25*(1.0+xi[0])
  sData.dhdxi[2,1] =  0.25*(1.0+xi[0])
  sData.dhdxi[3,1] =  0.25*(1.0-xi[0])

  sData.dhdxi[0,2] = -0.25*(1.0-xi[0])
  sData.dhdxi[1,2] = -0.25*(1.0+xi[0])
  sData.dhdxi[2,2] =  0.25*(1.0+xi[0])
  sData.dhdxi[3,2] =  0.25*(1.0-xi[0])

  return sData

#----------------------------------------------------------------------

def getShapePenta6 ( xi ):

  #Check the dimension of physical space
  if len(xi) != 3:
    raise NotImplementedError('3D only')

  #Initialise tuples
  
  sData = shapeData()
   
  sData.h     = empty( 6 )
  sData.dhdxi = empty( shape=(6,3) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = 0.25*(1.0-xi[0])*(1.0-xi[1])
  sData.h[1] = 0.25*(1.0+xi[0])*(1.0-xi[1])
  sData.h[2] = 0.25*(1.0+xi[0])*(1.0+xi[1])
  sData.h[3] = 0.25*(1.0-xi[0])*(1.0+xi[1])

  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.25*(1.0-xi[1])
  sData.dhdxi[1,0] =  0.25*(1.0-xi[1])
  sData.dhdxi[2,0] =  0.25*(1.0+xi[1])
  sData.dhdxi[3,0] = -0.25*(1.0+xi[1])

  sData.dhdxi[0,1] = -0.25*(1.0-xi[0])
  sData.dhdxi[1,1] = -0.25*(1.0+xi[0])
  sData.dhdxi[2,1] =  0.25*(1.0+xi[0])
  sData.dhdxi[3,1] =  0.25*(1.0-xi[0])

  sData.dhdxi[0,2] = -0.25*(1.0-xi[0])
  sData.dhdxi[1,2] = -0.25*(1.0+xi[0])
  sData.dhdxi[2,2] =  0.25*(1.0+xi[0])
  sData.dhdxi[3,2] =  0.25*(1.0-xi[0])

  return sData

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def getShapeHexa8 ( xi ):

  if len(xi) != 3:
    raise NotImplementedError('The isoparamatric coordinate should be 3D.')

  sData = shapeData()
  
  sData.h     = empty( 8 )
  sData.dhdxi = empty( shape=(8,3) )
  sData.xi    = xi

  #Calculate shape functions
  sData.h[0] = 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0-xi[2])
  sData.h[1] = 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[2])
  sData.h[2] = 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[2])
  sData.h[3] = 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[2])
  sData.h[4] = 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[2])
  sData.h[5] = 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0+xi[2])
  sData.h[6] = 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0+xi[2])
  sData.h[7] = 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[2])
 
  #Calculate derivatives of shape functions
  sData.dhdxi[0,0] = -0.125*(1.0-xi[1])*(1.0-xi[2])
  sData.dhdxi[1,0] =  0.125*(1.0-xi[1])*(1.0-xi[2])
  sData.dhdxi[2,0] =  0.125*(1.0+xi[1])*(1.0-xi[2])
  sData.dhdxi[3,0] = -0.125*(1.0+xi[1])*(1.0-xi[2])
  sData.dhdxi[4,0] = -0.125*(1.0-xi[1])*(1.0+xi[2])
  sData.dhdxi[5,0] =  0.125*(1.0-xi[1])*(1.0+xi[2])
  sData.dhdxi[6,0] =  0.125*(1.0+xi[1])*(1.0+xi[2])
  sData.dhdxi[7,0] = -0.125*(1.0+xi[1])*(1.0+xi[2])

  sData.dhdxi[0,1] = -0.125*(1.0-xi[0])*(1.0-xi[2])
  sData.dhdxi[1,1] = -0.125*(1.0+xi[0])*(1.0-xi[2])
  sData.dhdxi[2,1] =  0.125*(1.0+xi[0])*(1.0-xi[2])
  sData.dhdxi[3,1] =  0.125*(1.0-xi[0])*(1.0-xi[2])
  sData.dhdxi[4,1] = -0.125*(1.0-xi[0])*(1.0+xi[2])
  sData.dhdxi[5,1] = -0.125*(1.0+xi[0])*(1.0+xi[2])
  sData.dhdxi[6,1] =  0.125*(1.0+xi[0])*(1.0+xi[2])
  sData.dhdxi[7,1] =  0.125*(1.0-xi[0])*(1.0+xi[2])

  sData.dhdxi[0,2] = -0.125*(1.0-xi[0])*(1.0-xi[1])
  sData.dhdxi[1,2] = -0.125*(1.0+xi[0])*(1.0-xi[1])
  sData.dhdxi[2,2] = -0.125*(1.0+xi[0])*(1.0+xi[1])
  sData.dhdxi[3,2] = -0.125*(1.0-xi[0])*(1.0+xi[1])
  sData.dhdxi[4,2] =  0.125*(1.0-xi[0])*(1.0-xi[1])
  sData.dhdxi[5,2] =  0.125*(1.0+xi[0])*(1.0-xi[1])
  sData.dhdxi[6,2] =  0.125*(1.0+xi[0])*(1.0+xi[1])
  sData.dhdxi[7,2] =  0.125*(1.0-xi[0])*(1.0+xi[1])

  return sData

#----------------------------------------------------------------------

def getElemType( elemCoords ):
  
  nNel = elemCoords.shape[0]
  rank = elemCoords.shape[1]
  
  if rank == 1:
    if nNel == 2:
      return "Line2"
    elif nNel == 3:
      return "Line3"
    else:
      raise NotImplementedError('No 1D element with '+str(nNel)+' nodes available')
  elif rank == 2:
    if nNel == 3:
      return "Tria3"
    elif nNel == 4:
      return "Quad4"
    elif nNel == 6:
      return "Tria6"
    elif nNel == 8:
      return "Quad8"
    elif nNel == 9:
      return "Quad9"              
    else:
      raise NotImplementedError('No 2D element with '+str(nNel)+' nodes available')
  elif rank == 3:
    if nNel == 4:
      return "Tetra4"
    elif nNel == 6:
      return "Penta6"
    elif nNel == 8:
      return "Hexa8"
    else:
      raise NotImplementedError('No 3D element with '+str(nNel)+' nodes available')
  else:
    raise NotImplementedError('Rank must be 1,2 or 3')

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def tria_scheme( order ):

  if order == 1:
    xi     = [[1.0/3.0,1.0/3.0]]
    weight = [ 0.5 ]
  elif order == 3:
    r1 = 1.0/6.0
    r2 = 2.0/3.0
    
    xi = [[r1,r1],[r2,r1],[r1,r2]]
    
    w1 = 1.0/6.0
  
    weight = [w1,w1,w1]
  elif order == 7:
    r1 = 0.5*0.1012865073235
    r2 = 0.5*0.7974269853531
    r4 = 0.5*0.4701420641051
    r6 = 0.0597158717898
    r7 = 1.0/3.0

    xi = [[r1,r1],[r2,r1],[r1,r2],[r4,r6],[r4,r4],[r6,r4],[r7,r7]]

    w1 = 0.1259391805448
    w4 = 0.1323941527885
    w7 = 0.225

    weight = [ w1,w1,w1,w4,w4,w4,w7 ]
  
  return xi,weight
          
#-----------------------------------------------------------------------

def getIntegrationPoints( elemType , order , scheme ):

  xi     = []
  weight = []
  
  if elemType[:-1] == "Line":
    if elemType == "Line2":
      stdOrder = 2
    elif elemType == "Line3":
      stdOrder = 3
    xi,weight = gauss_scheme( stdOrder + order )
    xi = [float(a.real) for a in xi]

  elif elemType[:-1] == "Tria":
    orderArray = [1,3,7]
    if elemType == "Tria3":
      stdOrder = 0
    elif elemType == "Tria6":
      stdOrder = 1  
    xi,weight = tria_scheme( orderArray[stdOrder + order] )

  elif elemType[:-1] == "Quad":  
    if elemType == "Quad4":
      stdOrder = 2
    elif elemType == "Quad8" or elemType == "Quad9":
      stdOrder = 3  
    stdOrder += order

    ip,w  = gauss_scheme( stdOrder )
    
    for i in range(stdOrder):
      for j in range(stdOrder):
        xi.    append( [float(ip[i].real),float(ip[j].real)] )
        weight.append( w[i]*w[j] )
        
  elif elemType[:-1] == "Hexa":  
    if elemType == "Hexa8":
      stdOrder = 2
 
    stdOrder += order

    ip,w  = gauss_scheme( stdOrder )
    
    for i in range(stdOrder):
      for j in range(stdOrder):
        for k in range(stdOrder):
          xi.    append( [float(ip[i].real),float(ip[j].real),float(ip[k].real)] )
          weight.append( w[i]*w[j]*w[k] )

  return xi , weight

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def calcWeight( jac ):

  n = jac.shape

  if n[0] == n[1]:
    return abs(det(jac))
  elif n[0] == 2 and n[1] == 1:
    return sqrt(sum(sum(jac*jac)))

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def getElemShapeData( elemCoords , order = 0 , method = 'Gauss' , elemType = 'Default' ):

  elemData = elemShapeData()
  
  if elemType == 'Default':  
    elemType = getElemType( elemCoords )
    
  (intCrds,intWghts) = getIntegrationPoints( elemType , order , method )
    
  for xi,weight in zip( intCrds , intWghts ):    
    try:
      sData = eval( 'getShape'+elemType+'(xi)' )
    except:
      raise NotImplementedError('Unknown type :'+elemType)

    jac = dot ( elemCoords.transpose() , sData.dhdxi )

    if jac.shape[0] == jac.shape[1]:
      sData.dhdx = dot ( sData.dhdxi , inv( jac ) )

    sData.weight = calcWeight( jac ) * weight
    sData.x      = dot(sData.h,elemCoords)

    elemData.sData.append(sData)

  return elemData
