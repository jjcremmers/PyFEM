############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  Copyright (C) 2011-2024. The code is written in 2011-2012 by            #
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
import vtk

def storeNodes( grid , globdat ):

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

def storeElements( grid , globdat , elementGroup = "All" ):

    rank = globdat.nodes.rank

    for element in globdat.elements.iterElementGroup( elementGroup ):
        el_nodes = globdat.nodes.getIndices(element.getNodes())
    
        insertElement( grid , el_nodes , rank , element.family )
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def storeDofField( grid , data , globdat , dofTypes , label ):

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

def storeDofFields( grid , data , globdat ):

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
       
def storeNodeField( grid , data , globdat , name ):

    d = vtk.vtkDoubleArray();
    d.SetName( name );
    d.SetNumberOfComponents(1);
            
    for i,l in enumerate(data):
        d.InsertComponent( i , 0 , l )
                                                                 
    grid.GetPointData().AddArray( d )
   
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
def storeElementField( grid , data , globdat , name ):
             
    d = vtk.vtkDoubleArray();
    d.SetName( label );
    d.SetNumberOfComponents(1);
            
    for i,l in enumerate(data):        
        d.InsertComponent( i , 0 , l )
                 
    grid.GetCellData().AddArray( d )        

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------  
     
def setCellNodes( cell , elemNodes ):

    '''
  
    '''
          
    for i,inod in enumerate(elemNodes):
        cell.GetPointIds().SetId(i,inod)
          
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
          
def insertElement( grid , elemNodes , rank , family ):

    '''
    Inserts an element 
    '''
  
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

def insert2Dcontinuum( grid , elemNodes ):

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
      
def insert3Dcontinuum( grid , elemNodes ):

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
    
def insert2Dinterface( grid , elemNodes ):
    
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

def insert3Dinterface( grid , elemNodes ):
    
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
       
def insert2Dsurface( grid , elemNodes ):        
        
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
        
def insert3Dsurface( grid , elemNodes ):        
        
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

def insertBeam( grid , elemNodes ):        
        
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

def insertShell( grid , elemNodes ):        
        
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
        
        
             
