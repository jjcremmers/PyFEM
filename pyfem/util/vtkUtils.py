############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by            #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since    #
#  then augmented and  maintained by Joris J.C. Remmers.                   #
#  All rights reserved.                                                    #
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

from numpy import zeros
from pyfem.util.BaseModule import BaseModule
from pyfem.util.dataStructures import GlobalData
import vtk
from typing import List, Union

def storeNodes(grid: vtk.vtkUnstructuredGrid, globdat: GlobalData) -> None:
    """
    Store node coordinates from the global data into a VTK grid.
    
    Converts node coordinates to 3D format (padding 2D coordinates with z=0.0)
    and inserts them into the VTK grid's point data.
    
    Args:
        grid: VTK grid object to store points in
        globdat: Global data object
    """

    points = vtk.vtkPoints()
        
    for nodeID in list(globdat.nodes.keys()):
        crd = globdat.nodes.getNodeCoords(nodeID)

        crd1 = zeros(3)            
            
        if len(crd) == 2:
            crd1[:2] = crd
            crd1[2]  = 0.0
        else:
            crd1 = crd
          
        points.InsertNextPoint(crd1)

    grid.SetPoints(points) 
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def storeElements(grid: vtk.vtkUnstructuredGrid, globdat: GlobalData, elementGroup: str = "All") -> None:
    """
    Store elements from the global data into a VTK grid.
    
    Iterates through elements in the specified group and inserts them into
    the VTK grid with appropriate cell types based on element family and rank.
    
    Args:
        grid: VTK grid object to store cells in
        globdat: Global data object containing element information
        elementGroup: Name of element group to store. Defaults to "All"
    """

    rank = globdat.nodes.rank

    for element in globdat.elements.iterElementGroup( elementGroup ):
        el_nodes = globdat.nodes.getIndices(element.getNodes())
    
        insertElement( grid , el_nodes , rank , element.family )
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def storeDofField(grid: vtk.vtkUnstructuredGrid, data, globdat: GlobalData, dofTypes: Union[List[str], str], label: str) -> None:
    """
    Store a degree-of-freedom field as point data in the VTK grid.
    
    Creates a VTK array containing values for specified DOF types at each node.
    Missing DOF types are filled with zeros.
    
    Args:
        grid: VTK grid object to add array to
        data: Array of DOF values indexed by global DOF numbers
        globdat: Global data object containing DOF information
        dofTypes: List of DOF type strings (e.g., ["u", "v", "w"]) or single string
        label: Name for the VTK array
    """

    d = vtk.vtkDoubleArray()
    d.SetName( label )
    d.SetNumberOfComponents(len(dofTypes))
          
    i = 0
              
    for nodeID in list(globdat.nodes.keys()):
        j = 0
        for dispDof in dofTypes:
            if dispDof in globdat.dofs.dofTypes:
                d.InsertComponent( i , j , data[globdat.dofs.getForType(nodeID,dispDof)] )
            else:
                d.InsertComponent( i , j , 0.0 )
            j+=1
        i+=1      
        
    grid.GetPointData().AddArray( d ) 
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def storeDofFields(grid: vtk.vtkUnstructuredGrid, data, globdat: GlobalData) -> None:
    """
    Store all available DOF fields into the VTK grid.
    
    Automatically detects and stores displacement, temperature, and phase fields
    if they exist in the global data.
    
    Args:
        grid: VTK grid object to add arrays to
        data: Array of all DOF values
        globdat: Global data object containing DOF information
    """

    dofTypes = [ [ "u", "v", "w" ] , "temp" , "phase" ]
    labels   = [ "displacements" , "temperature" , "phase" ]
    
    for dofs,name in zip(dofTypes,labels):
        
        if type(dofs) == list:
            checkDof = dofs[0]
        else:
            checkDof = dofs
            
        if checkDof in globdat.dofs.dofTypes: 
            storeDofField( grid , data , globdat , dofs , name )
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
       
def storeNodeField(grid: vtk.vtkUnstructuredGrid, data, globdat: GlobalData, name: str) -> None:
    """
    Store a scalar field defined at nodes as point data in the VTK grid.
    
    Args:
        grid: VTK grid object to add array to
        data: Array of scalar values, one per node
        globdat: Global data object (for consistency with other functions)
        name: Name for the VTK array
    """

    d = vtk.vtkDoubleArray();
    d.SetName( name );
    d.SetNumberOfComponents(1);
            
    for i,l in enumerate(data):
        d.InsertComponent( i , 0 , l )
                                                                 
    grid.GetPointData().AddArray( d )
   
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
def storeElementField(grid: vtk.vtkUnstructuredGrid, data, globdat: GlobalData, name: str) -> None:
    """
    Store a scalar field defined at elements as cell data in the VTK grid.
    
    Args:
        grid: VTK grid object to add array to
        data: Array of scalar values, one per element
        globdat: Global data object (for consistency with other functions)
        name: Name for the VTK array
    """

    d = vtk.vtkDoubleArray();
    d.SetName( name );
    d.SetNumberOfComponents(1);
            
    for i,l in enumerate(data):        
        d.InsertComponent( i , 0 , l )
                 
    grid.GetCellData().AddArray( d )        

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------  
     
def setCellNodes(cell: vtk.vtkCell, elemNodes: List[int]) -> None:
    """
    Set the node IDs for a VTK cell.
    
    Maps the provided element node indices to the VTK cell's point IDs.
    
    Args:
        cell: VTK cell object to set nodes for
        elemNodes: List of node indices for the cell
    """

    '''
  
    '''
          
    for i,inod in enumerate(elemNodes):
        cell.GetPointIds().SetId(i,inod)
          
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
          
def insertElement(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int], rank: int, family: str) -> None:
    """
    Insert an element into the VTK grid.
    
    Routes element insertion to appropriate function based on element family
    (CONTINUUM, INTERFACE, SURFACE, BEAM, SHELL) and spatial rank (2D or 3D).
    
    Args:
        grid: VTK grid object to insert element into
        elemNodes: List of node indices for the element
        rank: Spatial dimension (2 or 3)
        family: Element family type
        
    Raises:
        NotImplementedError: If element family or rank combination is not supported
    """

    nNod = len(elemNodes)
  
    if family == "CONTINUUM":
        if rank == 2:
            insert2Dcontinuum( grid , elemNodes )               
        elif rank == 3:
            insert3Dcontinuum( grid , elemNodes )    
        else:
            raise NotImplementedError('Only 2D and 3D continuum elements.')
    elif family == "INTERFACE":
        if rank == 2:
            insert2Dinterface( grid , elemNodes )         
        elif rank == 3:
            insert3Dinterface( grid , elemNodes )        
        else:
            raise NotImplementedError('Only 2D and 3D interface elements.')
    elif family == "SURFACE":
        if rank == 2:
            insert2Dsurface( grid , elemNodes )                 
        elif rank == 3:
            insert3Dsurface( grid , elemNodes )                
    elif family == "BEAM":
        insertBeam( grid , elemNodes )        	     
    elif family == "SHELL":
        insertShell( grid , elemNodes )        
    else:
        raise NotImplementedError('Family of elements is not known.')
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def insert2Dcontinuum(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a 2D continuum element into the VTK grid.
    
    Supports 2, 3, 4, 6, 8, and 9 node elements. Higher-order elements are
    converted to linear elements using only corner nodes.
    
    Args:
        grid: VTK grid object to insert element into
        elemNodes: List of node indices for the element
        
    Raises:
        NotImplementedError: If element has unsupported number of nodes
    """

    nNod = len(elemNodes)
     
    if nNod == 2:
        cell = vtk.vtkLine()    
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() ) 
    elif nNod == 3:
        cell = vtk.vtkTriangle()    
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )     
    elif nNod == 4:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )       
    elif nNod == 6:
        cell = vtk.vtkTriangle()      
        setCellNodes( cell , elemNodes[0:6:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() ) 
    elif nNod == 8 or nNod == 9:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes[0:8:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
    else:
        print(nNod)
        raise NotImplementedError('Only 2, 3, 4, 6, 8, 9 node continuum elements in 2D.') 

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
def insert3Dcontinuum(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a 3D continuum element into the VTK grid.
    
    Supports 4, 5, 6, 8, and 16 node elements (tetrahedra, pyramids, wedges,
    and hexahedra). Higher-order elements are converted to linear elements.
    
    Args:
        grid: VTK grid object to insert element into
        elemNodes: List of node indices for the element
        
    Raises:
        NotImplementedError: If element has unsupported number of nodes
    """

    nNod = len(elemNodes)
    
    if nNod == 4:
        cell = vtk.vtkTetra()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )             
    elif nNod == 5:
        cell = vtk.vtkPyramid()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )            
    elif nNod == 6:
        cell = vtk.vtkWedge()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )             
    elif nNod == 8:
        cell = vtk.vtkHexahedron()
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )               
    elif nNod == 16:
        cell = vtk.vtkHexahedron()
        setCellNodes( cell , numpy.concatenate(elemNodes[0:8:2],elemNodes[8:16:2] ) ) 
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                         
    else:
        raise NotImplementedError('Only 4, 5, 6, 8 and 16 node continuum elements in 3D.')        
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
def insert2Dinterface(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a 2D interface element into the VTK grid.
    
    Represents interface elements as two separate line segments.
    
    Args:
        grid: VTK grid object to insert elements into
        elemNodes: List of 4 node indices (2 nodes per side)
        
    Raises:
        NotImplementedError: If element does not have 4 nodes
    """

    nNod = len(elemNodes)
                    
    if nNod == 4:
        cell = vtk.vtkLine() 
        setCellNodes( cell , elemNodes[0:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[2:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
    else:
        raise NotImplementedError('Only 4 node interface elements in 2D.')  
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def insert3Dinterface(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a 3D interface element into the VTK grid.
    
    Represents interface elements as two separate surfaces (triangles or quads).
    
    Args:
        grid: VTK grid object to insert elements into
        elemNodes: List of 6 (triangular) or 8 (quadrilateral) node indices
        
    Raises:
        NotImplementedError: If element does not have 6 or 8 nodes
    """

    nNod = len(elemNodes)
                    
    if nNod == 6:
        cell = vtk.vtkTria() 
        setCellNodes( cell , elemNodes[0:3] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[3:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )              
    elif nNod == 8:
        cell = vtk.vtkQuad() 
        setCellNodes( cell , elemNodes[0:4] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[4:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )      
    else:
        raise NotImplementedError('Only 6 and 8 node interface elements in 3D.')          
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
       
def insert2Dsurface(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a 2D surface element into the VTK grid.
    
    Args:
        grid: VTK grid object to insert element into
        elemNodes: List of 2 node indices
        
    Raises:
        NotImplementedError: If element does not have 2 nodes
    """

    nNod = len(elemNodes)
                    
    if nNod == 2:         
        cell = vtk.vtkLine() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
    else:
        raise NotImplementedError('Only 2 node surface elements in 2D.')       
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
        
def insert3Dsurface(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a 3D surface element into the VTK grid.
    
    Supports triangular and quadrilateral surface elements.
    
    Args:
        grid: VTK grid object to insert element into
        elemNodes: List of 3 (triangular) or 4 (quadrilateral) node indices
        
    Raises:
        NotImplementedError: If element does not have 3 or 4 nodes
    """

    nNod = len(elemNodes)   
    
    if nNod == 3:
        cell = vtk.vtkTriangle() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )              
    elif nNod == 4:
        cell = vtk.vtkQuad() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )            
    else:
        raise NotImplementedError('Only 3 and 4 node surface elements in 3D.')
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def insertBeam(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a beam element into the VTK grid.
    
    Supports 2-node (linear) and 3-node (quadratic) beam elements. For 3-node
    elements, only the end nodes are used.
    
    Args:
        grid: VTK grid object to insert element into
        elemNodes: List of 2 or 3 node indices
        
    Raises:
        NotImplementedError: If element does not have 2 or 3 nodes
    """

    nNod = len(elemNodes) 
    
    if nNod == 2:
        cell = vtk.vtkLine() 
        setCellNodes( cell , elemNodes )
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
    elif nNod == 3:
        cell = vtk.vtkLine()
        cell.GetPointIds().SetId(0,elemNodes[0])
        cell.GetPointIds().SetId(1,elemNodes[2])
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
    else:
        raise NotImplementedError('Only 2 and 3 node beam elements.')
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def insertShell(grid: vtk.vtkUnstructuredGrid, elemNodes: List[int]) -> None:
    """
    Insert a shell element into the VTK grid.
    
    Supports triangular and quadrilateral shell elements.
    
    Args:
        grid: VTK grid object to insert element into
        elemNodes: List of 3 (triangular) or 4 (quadrilateral) node indices
        
    Raises:
        NotImplementedError: If element does not have 3 or 4 nodes
    """

    nNod = len(elemNodes)         
                               
    if nNod == 3:
        cell = vtk.vtkTriangle() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
    elif nNod == 4:
        cell = vtk.vtkQuad() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                     
    else:
        raise NotImplementedError('Only 3 and 4 node shell elements.')



