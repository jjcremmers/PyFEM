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
import re
from pyfem.util.itemList import itemList

class ElementSet( itemList ):

  def __init__ ( self, nodes, props ):

    itemList.__init__( self )
    
    self.nodes  = nodes
    self.props  = props
    self.groups = {}

    ###############################
    # Create the material manager #
    ###############################

    #from pyfem.materials.MaterialManager import MaterialManager

    #if hasattr( props, 'material' ):
    #  self.matman = MaterialManager( props.materials )

  def __iter__ ( self ):

    elements = []

    for groupName in self.iterGroupNames():
      for element in self.iterElementGroup( groupName ):
        elements.append( element )
       
    return iter( elements )

  def getDofTypes ( self ):

    dofTypes = []

    for element in self:
      for dofType in element.dofTypes:
	      if dofType not in dofTypes:
	        dofTypes.append( dofType )

    return dofTypes
    
  def readFromFile( self, fname ):
    
    print("  Reading elements .............")

    fin = open( fname )
  
    while True:
      line = fin.readline()  
  
      if line.startswith('<Elements>') == True:
        while True:
          line = fin.readline()  

          if line.startswith('</Elements>') == True:
            return
        
          line = re.sub('\s{2,}',' ',line)
          a = line.split(';')
     
          for a0 in a[:-1]:
            b = a0.strip().split(' ')
            
            if b[0].startswith("//") or b[0].startswith("#"):
              break
            if len(b) > 1 and type(eval(b[0])) == int:
              self.add( eval(b[0]), eval(b[1]) , [eval(nodeID) for nodeID in b[2:]] )  

      elif line.startswith('gmsh') == True:
        ln = line.replace('\n','').replace('\t','').replace(' ','').replace('\r','').replace(';','')
        ln = ln.split('=',1)
        self.readGmshFile( ln[1][1:-1] )
        return

#
#
#

  def readGmshFile( self , fname ):

    fin = open( fname )

    while True:
    
      line = fin.readline()  
  
      if line.startswith('$MeshFormat') == True:
      
        while True:
          line = fin.readline()  
          line = re.sub('\s{2,}',' ',line)

          a = line.split(';')
          b = a[0].strip().split(' ')
      
          if eval(b[0]) < 2.0:
            print("error")
            sys.exit()
          
          break

      if line.startswith('$Elements') == True:
  
        nElements = eval(fin.readline())
    
        for i in range(nElements):
          line = fin.readline()
          line = re.sub('\s{2,}',' ',line)
          b    = line.strip().split(' ')

#          if len(b) > 1 and type(eval(b[0])) == int:
          if len(b) == 8 and type(eval(b[0])) == int:        
            self.add( eval(b[0]), "ContElem" , [eval(nodeID) for nodeID in b[5:]] )
      
      if line.startswith('$EndElements'):
        return 

#
#  add element
#
  def add ( self, ID, modelName, elementNodes ):  

    #Check if the model exists
    if not hasattr( self.props, modelName ):
      RuntimeError('Missing properties for model ' + modelName)

    modelProps = getattr( self.props, modelName )

    #Check if the model has a type
    if not hasattr( modelProps, 'type' ):
      RuntimeError('Missing type for model ' + modelName)
      
    modelType = getattr( modelProps, 'type' )

    #Load the element (Python 2.x)
    #cmdStr = 'from pyfem.elements.' + modelType + ' import ' + modelType + ' as element' 
    #exec(cmdStr)
 
    modelProps.rank = self.nodes.rank

    element = getattr(__import__('pyfem.elements.'+modelType , globals(), locals(), modelType , 0 ), modelType )

    #Create the element
 
    elem = element( elementNodes , modelProps )

    #  Check if the node IDs are valid:

    for nodeID in elem.getNodes():
      if not nodeID in self.nodes:
        raise RuntimeError('Node ID ' + str(nodeID) + ' does not exist')

    #  Add the element to the element set:

    itemList.add( self, ID, elem )

    #  Add the element to the correct group:

    self.addToGroup( modelName, ID )

  def addToGroup( self, modelType, ID ):

    if modelType not in self.groups:
      self.groups[modelType] = [ID]
    else:
      self.groups[modelType].append( ID )

  def addGroup ( self, groupName,  groupIDs ):
    self.groups[groupName] = groupIDs

  def iterGroupNames ( self ):
    return self.groups

  def iterElementGroup ( self, groupName ):
    if groupName == "All":
      return iter( self )
    else:
      return iter( self.get( self.groups[groupName] ) )

  def elementGroupCount( self, groupName ):
    if groupName == "All":
      return len(self)
    else:
      return len(self.groups[groupName])

  def commitHistory ( self ):

    for element in list(self.values()):
      element.commitHistory()
