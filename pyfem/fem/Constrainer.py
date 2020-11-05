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

class Constrainer:

  '''Constrainer class'''
  
  def __init__( self , nDofs , name = "Main" ):
  
    self.nDofs           = nDofs
    self.constrainData   = {}
    self.name            = name
    
    self.constrainedDofs = {}
    self.constrainedVals = {}
    self.constrainedFac  = {}
    
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
           
  def addConstraint(self,dofID,val,label):
  
    self.constrainData[dofID] = val   
    
    if dofID in self.constrainedDofs[label]:
      self.setFactorForDof( val , dofID, label )
      return
           
    self.constrainedDofs[label].append( dofID )
    
    if type(val) is list:
      if len(val) == 3:
        self.constrainedVals[label].append( val[0] )
    else:
      self.constrainedVals[label].append( val )
            
    if type(val) is list:
      if len(val) == 3:
        self.constrainData[val[1]] = "S"
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def flush ( self ):
  
    '''Returns the constraints matrix using the class member arrays'''
    
    row = []
    col = [] 
    val = [] 

    j = 0
    
    for i in range(self.nDofs):
      
      if i in self.constrainData:
        if type(self.constrainData[i]) is not list:
          #this a pure slave            
          continue
        else:
          slaves = self.constrainData[i]
                              
          row.append(slaves[1])
          col.append(j)
          val.append(slaves[2])
                    
      row.append(i)
      col.append(j)
      val.append(1.)
        
      j+=1
        
    self.C = coo_matrix((val,(row,col)), shape=(self.nDofs,j))


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def addConstrainedValues( self , a ):
  
    for name in self.constrainedDofs.keys():            
      a[self.constrainedDofs[name]] += self.constrainedFac[name] * array(self.constrainedVals[name])
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def setConstrainedValues( self , a ):
  
    for name in self.constrainedDofs.keys():             
      a[self.constrainedDofs[name]] = self.constrainedFac[name] * array(self.constrainedVals[name])
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def setConstrainFactor( self , fac , loadCase = "All_" ):

    if loadCase == "All_":
      for name in self.constrainedFac.keys():
        self.constrainedFac[name] = fac
    else:
      self.constrainedFac[loadCase] = fac
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------      

  def setFactorForDof( self , fac , dofID , label ):
    
    idx = self.constrainedDofs[label].index(dofID)
    self.constrainedVals[label][idx] = fac
    
#---------------------------------------------------
#
#--------------------------------------------------------

  def slaveCount( self ):
  
    counter = 0
    for name in self.constrainedFac.keys():
      counter += len(self.constrainedDofs[name])
      
    return counter  
      
