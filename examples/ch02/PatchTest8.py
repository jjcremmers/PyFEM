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
############################################################################
#  Description: The Python file of the example presented in section        #
#               2.6 of the book, pages 53--62. In this file, the mesh is   #
#               constructed using 8 node serendipity elements.             #
#                                                                          #
#  Use:         python PatchTest.py                                        #
############################################################################

from numpy import zeros,dot,ix_,array
from pyfem.util.shapeFunctions  import getElemShapeData
import scipy.linalg

#------------

def getDofs( nodes ):

  n = 2*len(nodes)

  dofs = zeros( n , dtype=int )

  dofs[0:n:2]=2*nodes
  dofs[1:n:2]=2*nodes+1

  return dofs

#--------------------

def getBmatrix( dhdx ):

  n = 2*len(dhdx)

  b = zeros( shape=( 3 , n ) )

  for i,dp in enumerate(dhdx):
    b[0,i*2  ] = dp[0]
    b[1,i*2+1] = dp[1]
    b[2,i*2  ] = dp[1]
    b[2,i*2+1] = dp[0]
  
  return b

#----------------------

nu = 0.25
E = 1.e6
D = zeros( shape=(3,3) )

D[0,0] = E/(1-nu*nu)
D[0,1] = D[0,0]*nu
D[1,0] = D[0,1]
D[1,1] = D[0,0]
D[2,2] = E/(2*(1+nu)	)

#Node coordinates

coords = zeros( shape=(20,2) )

coords[ 0,:] = [0.0 ,0.0 ]
coords[ 1,:] = [0.12,0.0 ]
coords[ 2,:] = [0.24,0.0 ]
coords[ 3,:] = [0.24,0.06]
coords[ 4,:] = [0.24,0.12]
coords[ 5,:] = [0.12,0.12]
coords[ 6,:] = [0.0 ,0.12]
coords[ 7,:] = [0.0 ,0.06]
coords[ 8,:] = [0.04,0.02]
coords[ 9,:] = [0.11,0.025]
coords[10,:] = [0.18,0.03]
coords[11,:] = [0.17,0.055]
coords[12,:] = [0.16,0.08]
coords[13,:] = [0.12,0.08]
coords[14,:] = [0.08,0.08]
coords[15,:] = [0.06,0.05]
coords[16,:] = [0.02,0.01]
coords[17,:] = [0.21,0.015]
coords[18,:] = [0.20,0.10]
coords[19,:] = [0.04,0.10]

#Elements

elems = zeros( shape=(5,8) , dtype=int )

elems[0,:] = [ 0 , 1 ,  2 , 17 , 10 ,  9 ,  8 , 16 ]
elems[1,:] = [ 2 , 3 ,  4 , 18 , 12 , 11 , 10 , 17 ]
elems[2,:] = [ 4 , 5 ,  6 , 19 , 14 , 13 , 12 , 18 ]
elems[3,:] = [ 6 , 7 ,  0 , 16 ,  8 , 15 , 14 , 19 ]
elems[4,:] = [ 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 ]

presNodes = array([0,1,2,3,4,5,6,7])

presInds  = getDofs( presNodes )
presVals  = zeros( len(presInds) )

upres = lambda crd : 1e-3*(crd[0]+crd[1]/2)
vpres = lambda crd : 1e-3*(crd[1]+crd[0]/2)

presVals[2*presNodes]   = [upres(crd) for crd in coords[presNodes,:] ]
presVals[2*presNodes+1] = [vpres(crd) for crd in coords[presNodes,:] ]

#----------------------------------------------------------------------
#  Calculate K
#----------------------------------------------------------------------

totDof = 2*len(coords)

K = zeros( shape=( totDof , totDof ) )

for elem in elems:
  elemDofs = getDofs(elem)

  sData = getElemShapeData( coords[elem,:] )

  for iData in sData:
    b = getBmatrix( iData.dhdx ) 
    K[ix_(elemDofs,elemDofs)] += dot ( b.transpose() , dot ( D , b ) ) * iData.weight

#----------------------------------------------------------------------
#  Solve Ka=f
#----------------------------------------------------------------------

consDof = len( presInds )
 
C = zeros( shape=(totDof,totDof-consDof) )

j = 0
    
for i in range(totDof):     
  if i in presInds:
    continue
  C[i,j] = 1.
  j+=1

a = zeros(totDof)

a[presInds] = presVals

Kcons = dot( dot( C.transpose(), K ), C )
fcons = dot( C.transpose(), dot( K , -a ) )

acons = scipy.linalg.solve( Kcons, fcons )

a = dot( C, acons )

a[presInds] = presVals

#----------------------------------------------------------------------
#  Calculate stresses and internal forces
#----------------------------------------------------------------------

fint        = zeros( totDof )
nodalStress = zeros( shape=(len(coords),3) )
nodalCount  = zeros( len(coords) )

for elem in elems:
  elemDofs = getDofs( elem )
  sData = getElemShapeData( coords[elem,:] )
 
  for iData in sData:
    b = getBmatrix( iData.dhdx )

    strain = dot( b , a[elemDofs] )
    stress = dot( D , strain )

    fint[elemDofs] += dot(b.transpose(),stress)*iData.weight
    
    nodalStress[elem,:] += stress
    nodalCount [elem]   += 1;

#----------------------------------------------------------------------
#  Print output
#----------------------------------------------------------------------

print(' Node|  d[x]       d[y]       |', end=' ')
print(' fint[x]    fint[y]    |', end=' ')
print(' sigma[xx]   sigma[yy]   tau[xy]')
print('---------------------------------------------------------------------------')

for i in range(len(coords)):
  print(' %3i | %10.3e %10.3e' % (i,a[2*i],a[2*i+1]), end=' ')
  print(' | %10.3e %10.3e' % (fint[2*i],fint[2*i+1]), end=' ')
  print(' | %10.3e' % (nodalStress[i,0]/nodalCount[i]), end=' ')
  print(' %10.3e' % (nodalStress[i,1]/nodalCount[i]), end=' ')
  print(' %10.3e' % (nodalStress[i,2]/nodalCount[i]))
