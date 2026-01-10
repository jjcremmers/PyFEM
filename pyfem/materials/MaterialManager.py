# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from importlib import import_module
from pyfem.util.dataStructures import Properties
from numpy import zeros
import copy

class MaterialManager ( list ):

    def __init__ ( self, matProps ):

        if hasattr( matProps,"type"):
            matType = matProps.type

            try:
                material = import_module(f"pyfem.materials.{matType}")
            except ModuleNotFoundError as e:
                raise ImportError(
                    f"Solver module 'pyfem.materials.{matType}' not found. "
                    f"Check the 'type' in your input file."
                ) from e    

            try:
                material_cls = getattr(material, matType)           
            except AttributeError as e:
                raise ImportError(
                    f"Class '{matType}' not found in module 'pyfem.materials.{matType}'. "
                    f"Ensure the class name matches the file name."
                ) from e

            self.material = material_cls(matProps)
            
            self.matlist     = []
            self.matProps    = matProps
            self.iSam        = -1
            self.failureFlag = False
    
        if hasattr(matProps,'failureType'):
    
            failureType = matProps.failureType
      
            try:
                failure = import_module(f"pyfem.materials.{failureType}")
            except ModuleNotFoundError as e:
                raise ImportError(
                    f"Solver module 'pyfem.materials.{failureType}' not found. "
                    f"Check the 'type' in your input file."
                ) from e    

            try:
                failure_cls = getattr(failure, failureType)           
            except AttributeError as e:
                raise ImportError(
                    f"Class '{failureType}' not found in module 'pyfem.material.{failureType}'. "
                    f"Ensure the class name matches the file name."
                ) from e

            self.failure = failure_cls(matProps)

            #failure = getattr(__import__('pyfem.materials.'+failureType , \
            #  globals(), locals(), failureType , 0 ), failureType )
 
            #self.failure = failure( matProps )
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
            self.matlist.append(self.material)
        
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