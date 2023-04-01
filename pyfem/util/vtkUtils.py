############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by            #
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

from pyfem.util.BaseModule import BaseModule
from pyfem.util.logger   import getLogger

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class vtkWriter():

  def __init__( self ):
  
  def storeNodes( self , nodes ):
  
  def storeElems( self , elems ):
  
  def addNodalData( self , data , label ):
  
  def addElementData( self , data , label ):
  
  def write( self , fileName ):
  
    vtkfile = open( fileName ,'w' )

    vtkfile.write('<?xml version="1.0"?>\n')
    vtkfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')

    vtkfile.write('<UnstructuredGrid>\n')

    vtkfile.write('<Piece NumberOfPoints="'+str(len(globdat.nodes))+'" NumberOfCells="')

    vtkfile.write(str(globdat.elements.elementGroupCount( self.elementGroup))+'">\n')

    vtkfile.write('<PointData>\n')

    vtkfile.write('<DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii" >\n')
  
class MeshWriter ( BaseModule ):
 
  def __init__( self , props , globdat ):
	
    self.prefix       = globdat.prefix
    self.elementGroup = "All"
    self.k            = 0
    self.interval     = 1
    self.extraFields  = []

    BaseModule.__init__( self , props )
    
    if type(self.extraFields) is str:
      self.extraFields = [self.extraFields]

  def run( self , props , globdat ):
    
    if not globdat.solverStatus.cycle%self.interval == 0:
      return

    logger.info("Writing mesh .................")

    dim = globdat.state.ndim    

    if dim == 1:
      self.writeCycle( globdat.state , props , globdat )
    elif dim == 2:
      for state in globdat.state.transpose():
        self.writeCycle( state , props , globdat )

    self.writePvd()

#
#
#

  def writeCycle( self , state , props , globdat ):

    vtkfile = open( self.prefix + '-' + str(self.k) + '.vtu' ,'w' )

    vtkfile.write('<?xml version="1.0"?>\n')
    vtkfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
    vtkfile.write('<UnstructuredGrid>\n')
    vtkfile.write('<Piece NumberOfPoints="'+str(len(globdat.nodes))+'" NumberOfCells="')
    vtkfile.write(str(globdat.elements.elementGroupCount( self.elementGroup))+'">\n')
    vtkfile.write('<PointData>\n')
    vtkfile.write('<DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii" >\n')
	
    dispDofs = ["u","v","w"]

    for nodeID in list(globdat.nodes.keys()):
      for dispDof in dispDofs:
        if dispDof in globdat.dofs.dofTypes:
          vtkfile.write(str(state[globdat.dofs.getForType(nodeID,dispDof)])+' ')
        else: 
          vtkfile.write(' 0.\n')
 
    vtkfile.write('</DataArray>\n')
    
    for field in self.extraFields:
      vtkfile.write('<DataArray type="Float64" Name="'+field+'" NumberOfComponents="1" format="ascii" >\n')
	
      for nodeID in list(globdat.nodes.keys()):      
        vtkfile.write(str(state[globdat.dofs.getForType(nodeID,field)])+' ')
  
      vtkfile.write('</DataArray>\n')
  
    for name in globdat.outputNames:
      stress = globdat.getData( name , list(range(len(globdat.nodes))) )

      vtkfile.write('<DataArray type="Float64" Name="'+name+'" NumberOfComponents="1" format="ascii" >\n')
      for i in range(len(globdat.nodes)):
        vtkfile.write( str(stress[i]) + " \n" )

      vtkfile.write('</DataArray>\n')
	
    vtkfile.write('</PointData>\n')
    vtkfile.write('<CellData>\n')
    vtkfile.write('</CellData>\n')
    vtkfile.write('<Points>\n')
    vtkfile.write('<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n')
  
    for nodeID in list(globdat.nodes.keys()):
      crd = globdat.nodes.getNodeCoords(nodeID)
      if len(crd) == 2:
        vtkfile.write( str(crd[0]) + ' ' + str(crd[1]) + " 0.0\n" )
      else:
        vtkfile.write( str(crd[0]) + ' ' + str(crd[1]) + ' ' + str(crd[2]) + "\n" )
    
    vtkfile.write('</DataArray>\n')
    vtkfile.write('</Points>\n')
    vtkfile.write('<Cells>\n')
    vtkfile.write('<DataArray type="Int64" Name="connectivity" format="ascii">\n')

    #--Store elements-----------------------------

    rank = globdat.nodes.rank

    for element in globdat.elements.iterElementGroup( self.elementGroup ):
      el_nodes = globdat.nodes.getIndices(element.getNodes())

      if rank == 2:
        if len(el_nodes) == 3 or len(el_nodes) == 4:
          for node in el_nodes:
            vtkfile.write(str(node)+' ')
        elif len(el_nodes) == 6 or len(el_nodes) == 8:
          for node in el_nodes[::2]:
            vtkfile.write(str(node)+' ')

      elif rank == 3:
        if len(el_nodes) == 8:
          for node in el_nodes:
            vtkfile.write(str(node)+' ')
 
      vtkfile.write('\n')
  
    vtkfile.write('</DataArray>\n')
    vtkfile.write('<DataArray type="Int64" Name="offsets" format="ascii">\n')

    for i,element in enumerate(globdat.elements.iterElementGroup( self.elementGroup )):
      nNel = len(globdat.nodes.getIndices(element.getNodes()))

      if rank == 2 and nNel == 8:
        nNel = 4
      vtkfile.write(str(nNel*(i+1))+'\n')

    vtkfile.write('</DataArray>\n')
    vtkfile.write('<DataArray type="UInt8" Name="types" format="ascii" RangeMin="9" RangeMax="9">\n')

    for element in globdat.elements.iterElementGroup( self.elementGroup ):
      nNel = len(globdat.nodes.getIndices(element.getNodes()))

      if rank == 2:
        if nNel == 3 or nNel ==6:
          vtkfile.write('5\n')
        else:
          vtkfile.write('9\n')
      else:
        if nNel == 8:
          vtkfile.write('12\n')

    vtkfile.write('</DataArray>\n')
    vtkfile.write('</Cells>\n')
    vtkfile.write('</Piece>\n')
    vtkfile.write('</UnstructuredGrid>\n')
    vtkfile.write('</VTKFile>\n') 
  
    self.k = self.k+1

#----------------------------------------------------------------------
#  writePvd
#----------------------------------------------------------------------

  def writePvd( self ):

    f = open( self.prefix + '.pvd' ,'w' )

    f.write("<VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>\n")
    f.write("<Collection>\n")
  
    for i in range(self.k):
      f.write("<DataSet file='"+self.prefix+'-'+str(i)+".vtu' groups='' part='0' timestep='"+str(i)+"'/>\n")
   
    f.write("</Collection>\n")
    f.write("</VTKFile>\n")

    f.close()
