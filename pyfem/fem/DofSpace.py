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

from numpy import array, dot, zeros
import scipy.linalg
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import eigsh
from pyfem.util.itemList import itemList
from pyfem.util.fileParser import readNodeTable

class DofSpace:

  def __init__ ( self, elements ):

    self.dofTypes = elements.getDofTypes()
    self.dofs     = array( list(range( len(elements.nodes) * len(self.dofTypes))) ).reshape( ( len(elements.nodes), len(self.dofTypes) ) )
    self.nodes    = elements.nodes

    #Create the ID map
    self.IDmap = itemList()
    for ind,ID in enumerate(elements.nodes):
      self.IDmap.add( ID, ind )

    self.constrainedDofs = []
    self.constrainedVals = []
    self.constrainedFac  = 1.

  def __str__ ( self ):
    return str(self.dofs)

  def __len__ ( self ):
    return len(self.dofs.flatten())

  def setConstrainFactor( self , fac ):
    self.constrainedFac = fac
    
  def readFromFile( self, fname ):
    
    print("  Reading constraints ..........\n")

    '''   
    fin = open( fname )

    while True:
      line = fin.readline()  
  
      if line.startswith('<NodeConstraints>') == True:
        while True:
          line = fin.readline()  

          if line.startswith('</NodeConstraints>') == True:
            return
        
          a = line.strip().split(';')
      
          if len(a) == 2:
            b = a[0].split('=')
        
            if len(b) == 2:
              c = b[0].split('[')
              
              dofType = c[0]
              nodeID  = eval(c[1].split(']')[0])
    '''

    nodeTable = readNodeTable( fname , "NodeConstraints" )

    for data in nodeTable[0].data:
      self.constrain( data[1] , data[0] , data[2] )
              
  def constrain ( self, nodeID, dofTypes , val = 0. ):

    if not nodeID in self.nodes:
      raise RuntimeError('Node ID ' + str(nodeID) + ' does not exist')

    ind = self.IDmap.get( nodeID )

    if isinstance( dofTypes, str ):
      dofTypes = [dofTypes]

    #Check if the dofTypes exist
    for dofType in dofTypes:
      if dofType not in self.dofTypes:
        raise RuntimeError('DOF type "' + dofType + '" does not exist')
      
    for dofType in dofTypes:    
      self.constrainedDofs.append( self.dofs[ind,self.dofTypes.index(dofType)] )
      self.constrainedVals.append( val )

  def getForType ( self, nodeIDs, dofType ):
    return self.dofs[self.IDmap.get( nodeIDs ),self.dofTypes.index(dofType)]

  def get ( self, nodeIDs ):
    return self.dofs[self.IDmap.get(nodeIDs)].flatten()

  def getConstraintsMatrix ( self ):

    n_constrained = len( self.constrainedDofs )
    n             = len( self )

    row = list(range(n))
    col = zeros( n , dtype=int )
    val = zeros(n)

    j = 0
    
    for i in range(n):
      
      if i in self.constrainedDofs:
        continue

      col[i]=j
      val[i]=1.
  
      j+=1
    
    return coo_matrix((val,(row,col)), shape=(n,n-n_constrained))

  def solve ( self, A, b ):

    if len(A.shape) == 2:
      C = self.getConstraintsMatrix()

      a = zeros(len(self))
      a[self.constrainedDofs] = self.constrainedFac * array(self.constrainedVals)

      A_constrained = C.transpose() * (A * C )
      b_constrained = C.transpose()* ( b - A * a )

      x_constrained = spsolve( A_constrained, b_constrained )

      x = C * x_constrained

      x[self.constrainedDofs] = self.constrainedFac * array(self.constrainedVals)
    
    elif len(A.shape) == 1:
      x = b / A

      x[self.constrainedDofs] = self.constrainedFac * array(self.constrainedVals)
   
    return x
    
  def eigensolve( self, A , B , count=5 ):

    C = self.getConstraintsMatrix()

    A_constrained = dot( dot( C.transpose(), A ), C )
    B_constrained = dot( dot( C.transpose(), B ), C )

    eigvals , eigvecs = eigsh( A_constrained, count , B_constrained , sigma = 0. , which = 'LM' )

    x = zeros(shape=(len(self),count))

    for i,psi in enumerate(eigvecs.transpose()):
      x[:,i] = C * psi
      
    return eigvals,x

  def norm ( self, r ):

    C = self.getConstraintsMatrix()
    
    return scipy.linalg.norm( C.transpose() * r )
