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

import numpy as np
from numpy import append,repeat

#===============================================================================
#
#===============================================================================


class Contact:

  def __init__ ( self , props ):
  
    self.flag = False
    self.dispDofs = ["u","v"]
    self.centre = [1.,1.]
    self.direction = [0.0,0.0]
    self.radius = 10.
    self.penalty = 1.e6
    
    if hasattr( props , 'contact' ):
      if hasattr( props.contact , 'type' ):
        self.type = props.contact.type
        self.flag = True
        self.centre    = np.array(props.contact.centre)
        self.radius    = props.contact.radius
        self.penalty   = props.contact.penalty
        self.direction = np.array(props.contact.direction)

#-------------------------------------------------------------------------------
#  checkContact   (with flag)
#-------------------------------------------------------------------------------
       
  def checkContact ( self , row , val , col , B , globdat ):
  
    if not self.flag:
      return
      
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
        
        B[idofs] += -self.penalty * overlap * normal
        
        mat = self.penalty * np.outer( normal , normal )

        row = append(row,repeat(idofs,len(idofs)))
        
        for i in range(len(idofs)):
          col=append(col,idofs)        

        val = append(val,mat.reshape(len(idofs)*len(idofs)))
              
    return row , val , col
