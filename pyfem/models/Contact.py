# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import numpy as np
from numpy import append,repeat
from pyfem.models.BaseModel import BaseModel

#===============================================================================
#
#===============================================================================


class Contact( BaseModel):

    def __init__ ( self , props , globdat ):
  
        BaseModel.__init__( self , props , globdat )
          
        self.flag = False
        self.dispDofs = ["u","v"]
        self.centre = [1.,1.]
        self.direction = [0.0,0.0]
        self.radius = 10.
        self.penalty = 1.e6
            
        if self.type == "sphere":
          self.dispDofs = ["u","v","w"] 
          
        self.flag      = True
        self.centre    = np.array(props.centre)
        self.radius    = props.radius
        self.penalty   = props.penalty
        self.direction = np.array(props.direction)

#-------------------------------------------------------------------------------
#  checkContact   (with flag)
#-------------------------------------------------------------------------------


    def run( self , props , globdat ):
          
        centre = self.centre + globdat.lam * self.direction

        for nodeID in list(globdat.nodes.keys()):
            crd = globdat.nodes.getNodeCoords(nodeID)
      
            idofs = globdat.dofs.getForTypes([nodeID],self.dispDofs)
   
            crd += globdat.state[idofs]
      
            ds = crd - centre
      
            dsnorm  = np.linalg.norm(ds)
            overlap = self.radius - dsnorm
                 
            if overlap > 0:
                normal = ds / dsnorm
        
                globdat.B[idofs] += -self.penalty * overlap * normal
        
                print("Contqact",crd,self.penalty,overlap)
                mat = self.penalty * np.outer( normal , normal )

                globdat.row = append(globdat.row,repeat(idofs,len(idofs)))
        
                for i in range(len(idofs)):
                    globdat.col=append(globdat.col,idofs)        

                globdat.val = append(globdat.val,mat.reshape(len(idofs)*len(idofs)))
