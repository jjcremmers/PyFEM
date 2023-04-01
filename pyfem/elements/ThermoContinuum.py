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

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from numpy import zeros, dot, outer, eye

class ThermoContinuum( Element ):
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.rank = props.rank

    self.dofTypes = [ 'temp' ]
    
    self.D     = self.material.heatConductivity*eye(self.rank)
    self.capac = self.material.heatCapacity
    
    if self.rank == 2:    
      self.labels = [ "qx" , "qy" ]
    elif self.rank == 3:
      self.labels = [ "qx" , "qy" , "qz" ]
    
    self.transient = True
    self.theta = 1.0
      
  def __type__ ( self ):
    return name
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):
       
    sData = getElemShapeData( elemdat.coords )
    
    nDof = len(elemdat.coords)
    
    temp0 = elemdat.state - elemdat.Dstate
    
    if self.transient:
      ctt      = zeros(shape=(nDof,nDof))
      invdtime = 1.0/self.solverStat.dtime
                       
    for iInt,iData in enumerate(sData):        
      gradTemp = dot( iData.dhdx.transpose() , elemdat.state )
            
      elemdat.stiff += \
        dot ( iData.dhdx , dot( self.D , iData.dhdx.transpose() ) ) * iData.weight
     
      if self.transient:
        ctt += self.capac * outer( iData.h , iData.h ) * iData.weight
              
      self.appendNodalOutput( self.labels , dot(self.D,gradTemp) ) 
    
    if self.transient:  
      ktt0 = invdtime * ctt - elemdat.stiff * ( 1.0-self.theta )
      
      elemdat.stiff *= self.theta
      
      elemdat.stiff += invdtime * ctt 
        
    elemdat.fint += dot ( elemdat.stiff , elemdat.state )
      
    if self.transient:
      elemdat.fint += -dot ( ktt0 , temp0 )
     
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
     
    nDof = self.dofCount()
    elemdat.stiff = zeros( shape=(nDof,nDof) )
    
    self.getTangentStiffness( elemdat )
