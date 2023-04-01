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

from pyfem.util.dataStructures import Properties
from numpy import zeros
import copy

class MaterialManager ( list ):

  def __init__ ( self, matProps ):

    if hasattr( matProps,"type"):
      matType = matProps.type

      self.material = getattr(__import__('pyfem.materials.'+matType , globals(), locals(), matType , 0 ), matType )
    
      self.matlist     = []
      self.matProps    = matProps
      self.iSam        = -1
      self.failureFlag = False
    
    if hasattr(matProps,'failureType'):
    
      failureType = matProps.failureType
      
      failure = getattr(__import__('pyfem.materials.'+failureType , \
        globals(), locals(), failureType , 0 ), failureType )
 
      self.failure = failure( matProps )
      self.failureFlag = True

  def reset( self ):

    self.iSam  = -1

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getStress ( self, kinematic , iSam = -1 ):

    '''
    
    '''
    
    if iSam == -1:
      self.iSam += 1
    else:
      self.iSam = iSam
            
    while self.iSam >= len(self.matlist):
      self.matlist.append(self.material( self.matProps ))
        
    self.mat = self.matlist[self.iSam]
    
    if self.mat.numericalTangent:
      self.mat.storeOutputFlag = True
      sigma1,tang = self.mat.getStress( kinematic)
                        
      self.mat.realNewHistory = copy.deepcopy(self.mat.newHistory)
      
      nStr = len(kinematic.strain)

      tan0 = zeros(shape=(nStr,nStr))
                   
      self.mat.storeOutputFlag = False

      for i in range(nStr):
        kin0 = copy.deepcopy(kinematic)
        kin0.strain[i]  += 1.0e-9
        kin0.dstrain[i] += 1.0e-9
        
        sigma,tang = self.mat.getStress( kin0)
        
        tan0[i,:] = (sigma - sigma1 )/ (kin0.strain[i]-kinematic.strain[i])
        
      result = (sigma1,tan0)
      
      self.mat.newHistory = copy.deepcopy(self.mat.realNewHistory)
      
    else:
      self.mat.storeOutputFlag = True
      result = self.mat.getStress( kinematic )
    
    if self.failureFlag:
      self.failure.check(result[0],kinematic)
      
    return result
    
  def getStressPiezo ( self, kinematic , elecField, iSam = -1 ):

    if iSam == -1:
      self.iSam += 1
    else:
      self.iSam = iSam
            
    while self.iSam >= len(self.matlist):
      self.matlist.append(self.material( self.matProps ))
        
    self.mat = self.matlist[self.iSam]
     
    result = self.mat.getStressPiezo( kinematic , elecField )
    
    if self.failureFlag:
      self.failure.check(result[0],kinematic, elecField)
      
    return result
    
  def outLabels( self ):
    return self.mat.outLabels

  def outData( self ):
    return self.mat.outData

  def getHistory( self , label ):
    return self.mat.getHistoryParameter( label )

  def commitHistory( self ):
  
    if hasattr(self,"matlist"):
      for mat in self.matlist:
        mat.commitHistory()
