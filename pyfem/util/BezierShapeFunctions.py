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

from numpy import array, dot, ndarray, empty, zeros , ones, real , sqrt
from scipy.linalg import norm , det , inv
from scipy.special.orthogonal import p_roots as gauss_scheme
import numpy

from .shapeFunctions import shapeData,elemShapeData,getIntegrationPoints,getElemType
            
#----------------------------------------------------------------------

def getBezierLine4 ( xi , C ):

  #Check the dimensions of the parametric space

  if type(xi) == 'numpy.float64':
    raise NotImplementedError('1D only')
  if C.shape[1] != 4:
    raise NotImplementedError('C needs to have 4 columns.')

  sData       = shapeData()
  
  #Set length of lists

#  sData.h     = empty( 4 )
#  sData.dhdxi = empty( shape=(1,4) )
  sData.xi    = xi
  
  B     = empty(4)
  dBdxi = empty( shape=(4,1))

  #Calculate shape functions

  B[0] = -0.125*(xi-1.)**3
  B[1] =  0.375*(xi-1.)**2*(xi+1.)
  B[2] = -0.375*(xi-1.)*(xi+1.)**2
  B[3] =  0.125*(xi+1.)**3

  #Calculate derivatives of shape functions

  dBdxi[0,0] = -0.375*(xi-1.)**2
  dBdxi[1,0] =  0.75*(xi-1.0)*(xi+1.0) + 0.375*(xi-1.)**2
  dBdxi[2,0] = -0.375*(1+xi)**2-0.75*(1+xi)*(xi-1)
  dBdxi[3,0] =  0.375*(xi+1.)**2

  sData.h     = dot( C , B )
  sData.dhdxi = dot( C , dBdxi )

  return sData

#-------------

def calcWeight( jac ):

  n = jac.shape

  if n[0] == n[1]:
    return det(jac)
  elif n[0] == 1 and n[1] == 2:
    return sqrt(sum(sum(jac*jac)))
          
#----------------------------------------------------------------------

def getElemBezierData( elemCoords , C , order=4 , method="Gauss" , elemType = 'default' ):

  elemData = elemShapeData()
    
  if elemType == 'default':
    elemType = getElemType( elemCoords )
    
  (intCrds,intWghts) = getIntegrationPoints( "Line3" , order , method )
    
  for xi,intWeight in zip( real(intCrds) , intWghts ):    
    try:
      sData = eval( 'getBezier'+elemType+'(xi,C)' )
    except:
      raise NotImplementedError('Unknown type :'+elemType)

    jac = dot ( sData.dhdxi.transpose() , elemCoords )

    if jac.shape[0] is jac.shape[1]:
      sData.dhdx   = (dot ( inv( jac ) , sData.dhdxi.transpose() )).transpose()
   
    sData.weight = calcWeight( jac ) * intWeight

    elemData.sData.append(sData)

  return elemData
