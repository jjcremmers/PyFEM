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

from pyfem.util.BaseModule import BaseModule
import h5py
import numpy as np

from pyfem.util.logger import getLogger

logger = getLogger()

class HDF5Writer( BaseModule ):

  def __init__ ( self, props , globdat ):

    self.prefix    = globdat.prefix
    self.extension = ".h5"
    
    self.dispDofs    = ["u","v","w"]
    self.extraFields = ["rx","ry","rz","temp","pres"]
    self.singleFile  = True
    
    BaseModule.__init__( self , props )
    
    if not hasattr( props , "interval" ):
      self.interval = 1
    
    if self.singleFile:  
      f = h5py.File(self.prefix+self.extension,"w")            
      f.attrs["cycleCount"] = 0
           
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def run( self , props , globdat ):
  
    cycle = globdat.solverStatus.cycle
    
    if cycle % self.interval == 0:
     
      logger.info("Writing hdf5 file ............")
      
      if self.singleFile:
        f = h5py.File(self.prefix+self.extension,"a")
            
        self.cycle = f.attrs["cycleCount"]
        
        self.cycle += 1
        
        f.attrs["cycleCount"] = self.cycle
        
        gName = "cycle"+str(self.cycle)
        f.create_group(gName)
        
        self.writeCycle(f[gName] , globdat )
        
      else:
        name = str(self.prefix + "_" + str(self.cycle) + self.extension)
      
        f = h5py.File(name, "w")
      
        f.attrs['fileFormat'] = "RNDF"
        f.attrs['version']    = 1.0
        
        self.writeCycle(f)
        
        self.cycle += 1
 
#----------------------------------------------------------------------------------
#  writeCycle
#----------------------------------------------------------------------------------
  
  def writeCycle( self , cdat , globdat ):  
      
    cdat.create_group("elements")
           
    elemCount    = []
    connectivity = []
      
    i0 = 0
    for elem in globdat.elements:    
      i0 = i0 + len(elem)
      elemCount.append(i0)
      connectivity.extend(elem)
       
    connectivity = np.array(globdat.nodes.getIndices(connectivity),dtype = int )
    elemCount    = np.array(elemCount,dtype = int )
    elemIDs      = np.array(globdat.elements.getIndices(),dtype = int )
    familyIDs    = np.array(globdat.elements.getFamilyIDs(),dtype = int )
              
    cdat["elements"].create_dataset("offsets",      elemCount.shape,    dtype='i', data=elemCount )
    cdat["elements"].create_dataset("connectivity", connectivity.shape, dtype='i', data=connectivity )
    cdat["elements"].create_dataset("elementIDs",   elemIDs.shape,      dtype='i', data=elemIDs )
    cdat["elements"].create_dataset("familyIDs",    familyIDs.shape,      dtype='i', data=familyIDs )    
    
    cdat.create_group("elementGroups")
      
    for key in globdat.elements.groups:
      elementIDs = np.array(globdat.elements.getIndices(globdat.elements.groups[key]), dtype = int )      
      cdat["elementGroups"].create_dataset( key , elementIDs.shape, dtype='i' , data=elementIDs )
            
    cdat.create_group("nodes")
            
    coordinates = []
      
    for nodeID in list(globdat.nodes.keys()):
      coordinates.append(globdat.nodes.getNodeCoords(nodeID))
         
    coordinates = np.array(coordinates,dtype = float)
    nodeIDs     = np.array(globdat.nodes.getIndices(),dtype = int )
              
    cdat["nodes"].create_dataset("coordinates", coordinates.shape, dtype='f', data=coordinates)   
    cdat["nodes"].create_dataset("nodeIDs",     nodeIDs.shape,     dtype='i', data=nodeIDs)       
              
    dofs = self.dispDofs[:coordinates.shape[1]]
          
    cdat.create_group("nodeGroups")      
       
    for key in globdat.nodes.groups:
      nodeIDs = np.array(globdat.nodes.getIndices(globdat.nodes.groups[key]), dtype = int )      
      cdat["nodeGroups"].create_dataset( key , nodeIDs.shape, dtype='i' , data=nodeIDs )
                
    cdat.create_group("nodeData")
           
    displacements = []

    for nodeID in list(globdat.nodes.keys()):
      d = []
      for dispDof in dofs:
        if dispDof in globdat.dofs.dofTypes:
          d.append(globdat.state[globdat.dofs.getForType(nodeID,dispDof)])
        else: 
          d.append(0.)
      displacements.append(d)
        
    displacements = np.array(displacements,dtype = float )
 
    cdat["nodeData"].create_dataset("displacements",displacements.shape, dtype='f', data=displacements)   
           
    for field in self.extraFields:     	              
      if field in globdat.dofs.dofTypes:  
        output = []
                     
        for nodeID in list(globdat.nodes.keys()):      
          output.append(globdat.state[globdat.dofs.getForType(nodeID,field)])
  
        output = np.array( output,dtype = float )
      
        cdat["nodeData"].create_dataset(field,output.shape, dtype='f', data=output)   
      
    for name in globdat.outputNames:      
      output = globdat.getData( name , list(range(len(globdat.nodes))) )
        
      output = np.array( output,dtype = float )
      
      cdat["nodeData"].create_dataset(name,output.shape, dtype='f', data=output)            
