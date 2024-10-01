################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2024. The code is written in 2011-2012 by                #
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
                
