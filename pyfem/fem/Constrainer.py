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

from numpy import array, dot, zeros
import scipy.linalg
from scipy.sparse import coo_matrix

from pyfem.util.logger   import getLogger

logger = getLogger()

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
  
    if dofID in self.constrainData:
      self.constrainData[dofID].append(val)
    else:
      self.constrainData[dofID] = [val]

    if (type(val) is list) and (len(val)==3):
      addVal = val[0]
    else:
      addVal = val

    if dofID in self.constrainedDofs[label]:
      self.setFactorForDof( addVal , dofID, label )     
      return
           
    self.constrainedDofs[label].append( dofID )
    
    self.constrainedVals[label].append( addVal )     

    logger.debug("TEST")

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def checkConstraints ( self , dofspace, nodeTables ):
    
    '''Checks tying relations between dofs'''
    for item in self.constrainData:

      dofInd = item

      for ilabel in self.constrainedDofs:
        if dofInd in self.constrainedDofs[ilabel]:
          label = ilabel

      removeItem = []

      for tiedItem in self.constrainData[item]:

        if (type(tiedItem) is list) and len(tiedItem)==3:
          Dat = tiedItem

          valSlave = Dat[0]
          masterDofID = Dat[1][0]
          facSlave = Dat[2]
          
          if masterDofID in self.constrainData:
            tempVal = []
            tempFac = []
            
            # Recursive loop until masterDofID not a list, but prescribed value         
            while masterDofID in self.constrainData:
              master = self.constrainData[masterDofID][0]
              if type(master) is list and len(master)==3:
                masterDofID = master[1][0]
                tempVal.append(master[0])
                tempFac.append(master[2])
              else:
                masterFin = master
                masterDofID = -1

            for iVal,iFac in reversed(list(zip(tempVal,tempFac))):
              masterFin += iVal + master*iFac
            
            self.addConstraint(dofInd,valSlave+masterFin*facSlave,label)
   
            removeItem.append(tiedItem)      

      for iRemove in removeItem:
        self.constrainData[item].remove(iRemove)  

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def flush ( self ):
  
    '''Returns the constraints matrix using the class member arrays'''
    
    row = []
    col = [] 
    val = [] 
    master = {}

    iCon = 0

    for iDof in range(self.nDofs):
      if iDof in self.constrainData:
        for item in self.constrainData[iDof]:
          if type(item) is not list:
            continue
          else:
            if item[1][0] in self.constrainData:
              #Something not checked correct in checkConstraint2
              raise RuntimeError('ERROR - Master of slave is a slave itself')
            else:
              master[iDof] = item
        
      else:          
        row.append(iDof)
        col.append(iCon)
        val.append(1.)

        iCon+=1
    
    # Assign correct slaves to masters
    for iSlave in master:
      fac = master[iSlave][2]
      masterDofID = master[iSlave][1][0]

      #Find column for free DOF of the Master
      rowID = row.index(masterDofID)
      col.append(rowID)
      row.append(iSlave)
      val.append(fac)
          
    self.C = coo_matrix((val,(row,col)), shape=(self.nDofs,iCon))

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
    self.constrainedVals[label][idx] += fac
    
#---------------------------------------------------
#
#--------------------------------------------------------

  def slaveCount( self ):
  
    counter = 0
    for name in self.constrainedFac.keys():
      counter += len(self.constrainedDofs[name])
      
    return counter  
      
