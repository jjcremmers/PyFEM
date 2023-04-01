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

import re
from pyfem.util.itemList import itemList
from pyfem.util.logger   import getLogger
from pyfem.util.dataStructures import solverStatus

logger = getLogger()

class ElementSet( itemList ):

  def __init__ ( self, nodes, props ):

    itemList.__init__( self )
    
    self.nodes  = nodes
    self.props  = props
    self.solverStat = solverStatus()
    self.groups = {}

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def __iter__ ( self ):

    elements = []

    for groupName in self.iterGroupNames():
      for element in self.iterElementGroup( groupName ):
        elements.append( element )
       
    return iter( elements )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def __repr__( self ):
    msg =  "Number of elements ......... %6d\n" % len(self)
    
    if len(self.groups) > 0:
      msg += "  Number of  groups .......... %6d\n" % len(self.groups)
      msg += "  -----------------------------------\n"
      msg += "    name                       #elems\n"
      msg += "    ---------------------------------\n"
      for name in self.groups:
        msg += "    %-16s           %6d\n" % (name,len(self.groups[name]))
    
    return msg

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getDofTypes ( self ):

    dofTypes = []

    for element in self:
      for dofType in element.dofTypes:
        if dofType not in dofTypes:
          dofTypes.append( dofType )

    return dofTypes
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def readFromFile( self, fname ):
    
    logger.info("Reading elements .............")
    
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
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def readGmshFile( self , fname ):

    import meshio
    
    mesh = meshio.read(fname,file_format="gmsh")
    
    elemID = 0

    for key in mesh.cell_sets_dict: 
      for typ in mesh.cell_sets_dict[key]:
        for idx in mesh.cell_sets_dict[key][typ]:
          iNodes = mesh.cells_dict[typ][idx]
          self.add( elemID , key , iNodes.tolist() )
          elemID = elemID + 1
                                 
#-------------------------------------------------------------------------------
#  add element
#-------------------------------------------------------------------------------

  def add ( self, ID, modelName, elementNodes ):  

    #Check if the model exists
    
    if hasattr( self.props, modelName ):
    
      modelProps = getattr( self.props, modelName )

      #Check if the model has a type
      if not hasattr( modelProps, 'type' ):
        raise RuntimeError('Missing type for model ' + modelName)
      
      modelType = getattr( modelProps, 'type' )
 
      modelProps.rank       = self.nodes.rank
      modelProps.solverStat = self.solverStat

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

  def addGroup ( self, groupName,  groupIDs ):
    self.groups[groupName] = groupIDs

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def iterGroupNames ( self ):
    return self.groups

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def iterElementGroup ( self, groupName ):
    if groupName == "All":
      return iter( self )
    elif isinstance(groupName, list):
      elems = []
      for name in groupName:
        elems += self.get(self.groups[name])
      return iter( elems )
    else:
      return iter( self.get( self.groups[groupName] ) )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def elementGroupCount( self, groupName ):
    if groupName == "All":
      return len(self)
    elif isinstance(groupName, list):
      length = 0;
      for name in groupName:
        length += len(self.groups[name])
      return length
    else:
      return len(self.groups[groupName])
      
#
#
#

  def getFamilyIDs(self):
  
    familyIDs = []
    fam = ["CONTINUUM","INTERFACE","SURFACE","BEAM","SHELL"]
    
    for elem in self:
      
      familyIDs.append(fam.index(elem.family))
      
    return familyIDs

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def commitHistory ( self ):

    for element in list(self.values()):
      element.commitHistory()
