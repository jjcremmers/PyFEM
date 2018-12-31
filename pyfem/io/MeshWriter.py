############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stabke version can be downloaded from the web-site:          #
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

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class MeshWriter ( BaseModule ):
 
  def __init__( self , props , globdat ):
	
    self.prefix       = globdat.prefix
    self.elementGroup = "All"
    self.k            = 0
    self.interval     = 1

    BaseModule.__init__( self , props )

  def run( self , props , globdat ):
    
    if not globdat.cycle%self.interval == 0:
      return

    print("  Writing mesh .................\n")

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
	
    for nodeID in list(globdat.nodes.keys()):
      for dofType in globdat.dofs.dofTypes:
        vtkfile.write(str(state[globdat.dofs.getForType(nodeID,dofType)])+' ')
   
      if len(globdat.dofs.dofTypes) == 2:
        vtkfile.write(' 0.\n')
 
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

    for element in globdat.elements.iterElementGroup( self.elementGroup ):
      el_nodes = globdat.nodes.getIndices(element.getNodes())

      if len(globdat.dofs.dofTypes) == 2:
        if len(el_nodes) == 3 or len(el_nodes) == 4:
          for node in el_nodes:
            vtkfile.write(str(node)+' ')
  
        elif len(el_nodes) == 6 or len(el_nodes) == 8:
          for node in el_nodes[::2]:
            vtkfile.write(str(node)+' ')

      if len(globdat.dofs.dofTypes) == 3:
        if len(el_nodes) == 8:
          for node in el_nodes:
            vtkfile.write(str(node)+' ')
  
      vtkfile.write('\n')
  
    vtkfile.write('</DataArray>\n')
    vtkfile.write('<DataArray type="Int64" Name="offsets" format="ascii">\n')

    for i,element in enumerate(globdat.elements.iterElementGroup( self.elementGroup )):
      el_nodes = globdat.nodes.getIndices(element.getNodes())
      vtkfile.write(str(len(el_nodes)*(i+1))+'\n')

    vtkfile.write('</DataArray>\n')
    vtkfile.write('<DataArray type="UInt8" Name="types" format="ascii" RangeMin="9" RangeMax="9">\n')

    for element in globdat.elements.iterElementGroup( self.elementGroup ):
      nNel = len(globdat.nodes.getIndices(element.getNodes()))
      if len(globdat.dofs.dofTypes) == 2:
        vtkfile.write('9\n')
      else:
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
