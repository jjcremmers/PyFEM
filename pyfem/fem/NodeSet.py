############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
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
from numpy import array
from pyfem.util.itemList import itemList
from pyfem.util.fileParser import getType
import re,sys,meshio

from pyfem.util.logger   import getLogger

logger = getLogger()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class NodeSet( itemList ):

  def __init__( self ):
    self.rank = -1
    self.groups = {}

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getNodeCoords( self, nodeIDs ):
    return array( self.get( nodeIDs ) )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def readFromFile( self, fname ):
    
    logger.info("Reading nodes ................")

    fin = open( fname , 'r' )
    
    line = fin.readline() 
    
    while line:    
      if line.replace(" ","").startswith('<Nodes>'):
        self.readNodalCoords( fin )

      if line.replace(" ","").startswith('gmsh'):
        ln = line.replace('\n','').replace('\t','').replace(' ','').replace('\r','').replace(';','')
        ln = ln.split('=',1)
        self.readGmshFile( ln[1][1:-1] )       
        break        
       
      line = fin.readline() 
      
    fin = open( fname , 'r' )
    
    line = fin.readline() 
    
    while line:  
      if line.replace(" ","").startswith('<NodeGroup'):
        if 'name' in line:
          label = line.split('=')[1].replace('\n','').replace('>','').replace(' ','').replace('\"','').replace('\'','')
          self.readNodegroup( fin , label )
       
      line = fin.readline()       
          
    for key in self.groups:
      self.groups[key] = list(set(self.groups[key]))
      
    print(self)
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def readGmshFile( self, fname ):
  
    mesh = meshio.read(fname,file_format="gmsh")

    self.rank = 2#len(mesh.points[0])
    
    for nodeID,p in enumerate(mesh.points):
      self.add(nodeID,p[:self.rank])

    for key in mesh.cell_sets_dict:
      for typ in mesh.cell_sets_dict[key]:
        for idx in mesh.cell_sets_dict[key][typ]:         
          iNodes = mesh.cells_dict[typ][idx]
          for nodeID in iNodes:
            self.addToGroup( key , nodeID )
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def addToGroup( self, modelType, ID ):

    if modelType not in self.groups:
      self.groups[modelType] = [ID]
    else:
      self.groups[modelType].append( ID )
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def __repr__( self ):
    msg =  "Nodeset contains %i nodes.\n" % len(self)
    
    if len(self.groups) > 0:
      msg += "-----------------------------------------\n"
      msg += "Number of nodegroups ... %i\n" % len(self.groups)
      
      for name in self.groups:
        msg += "  %s contains %i nodes \n" % (name,len(self.groups[name]))
    
    return msg
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def readNodalCoords( self , fin ):
    
    while True:
      line = fin.readline()  

      if line.replace(" ", "").startswith('</Nodes>'):
        return
  
      line = re.sub('\s{2,}',' ',line)
      a = line.split(';')
     
      for a in a[:-1]:
        b = a.strip().split(' ')
            
        if b[0].startswith("//") or b[0].startswith("#"):
          break
        if len(b) > 1 and type(eval(b[0])) == int:
          if self.rank == -1:
            self.rank = len(b)-1

          self.add( eval(b[0]), [eval(crd) for crd in b[1:]] ) 
          
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def readNodegroup( self , fin , key ):
    
    while True:
      line = fin.readline()
      
      if line.replace(" " ,"").startswith('</NodeGro'):
        return
        
      a = line.split()
      
      for b in a:
        if getType(b) == int:
          self.addToGroup( key ,b )
