# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util.BaseModule import BaseModule
import vtk

from pyfem.util.vtkUtils import ( insertElement,storeNodes,storeElements,
                                  storeDofFields,storeDofField,storeNodeField,
                                  storeElementField )
from numpy import zeros
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class MeshWriter ( BaseModule ):
 
    def __init__( self , props , globdat ):
	
        self.prefix       = globdat.prefix
        self.elementGroup = "All"
        self.interval     = 1
        self.extraFields  = []
        self.beam         = False
        self.interface    = False
        self.format       = "binary"

        BaseModule.__init__( self , props )
    
        if type(self.extraFields) is str:
            self.extraFields = [self.extraFields]
            
        self.vtufiles = []
        self.cycles   = []

    def run( self , props , globdat ):
    
        if not globdat.solverStatus.cycle%self.interval == 0:
            return
      
        self.writeHeader( globdat.solverStatus.cycle )

        dim = globdat.state.ndim    

        if dim == 1:
            self.writeCycle( globdat.state , props , globdat )
        elif dim == 2:
            for state in globdat.state.transpose():
                self.writeCycle( state , props , globdat )

        self.writePvd()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def writeCycle( self , state , props , globdat ):
  
        dispDofs = ["u","v","w"]
        cycle    = globdat.solverStatus.cycle
            
        writer = vtk.vtkXMLUnstructuredGridWriter()
  
        vtufile = self.prefix+'_t'+str(cycle)+".vtu"
            
        writer.SetFileName(vtufile)
      
        self.vtufiles.append(vtufile)
        self.cycles.  append(cycle)

        grid = vtk.vtkUnstructuredGrid()
        
        storeNodes( grid , globdat )
            
        #--Store elements-----------------------------
        
        storeElements( grid , globdat , self.elementGroup )
                      
        # -- Write nodedata
        
        storeDofFields( grid , state , globdat )
        
        # -- Write modes
                
        if hasattr( globdat , "eigenvecs" ):
            for iMod,eigenvecs in enumerate(globdat.eigenvecs.T):
                storeDofField( grid , eigenvecs , globdat , [ "u", "v", "w" ] , "mode"+str(iMod))        
                  
        # ------
               
        for name in globdat.outputNames:
            data = globdat.getData( name , list(range(len(globdat.nodes))) )
            
            storeNodeField( grid , data , globdat , name )
        
        # -- Write elemdata

        if hasattr( globdat , "elementData" ): 
            elemData = globdat.elementData
                        
            for label in elemData.outputNames:
                data    = getattr( elemData , label )
                         
                storeElementField( grid , data , globdat , label )
            
        
        writer.SetInputData(grid)
        
        if self.format == 'binary':
            writer.SetDataModeToBinary()
        else:
            writer.SetDataModeToAscii()

        writer.Write()    

#-------------------------------------------------------------------------------
#  writePvd
#-------------------------------------------------------------------------------

    def writePvd( self ):

        f = open( self.prefix + '.pvd' ,'w' )

        f.write("<VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>\n")
        f.write("<Collection>\n")
  
        for cycle,fileName in zip(self.cycles,self.vtufiles):
            f.write("<DataSet file='"+fileName+"' groups='' part='0' timestep='"+str(cycle)+"'/>\n")
   
        f.write("</Collection>\n")
        f.write("</VTKFile>\n")

        f.close()
                
