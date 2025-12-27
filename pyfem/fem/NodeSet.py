################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
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

from numpy import array
from typing import Any, Dict, Iterable, Iterator, List, TextIO, Union
from pyfem.util.itemList import itemList
from pyfem.util.fileParser import getType
import re, sys

from pyfem.util.logger import getLogger

logger = getLogger()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class NodeSet(itemList):
    """Container for nodes and node groups.

    Maintains nodal coordinates, dimensionality (`rank`), and groupings read
    from legacy `.pro` files or Gmsh meshes. Provides convenience methods for
    retrieving coordinates and iterating group contents.
    """

    def __init__(self) -> None:
        self.rank: int = -1
        self.groups: Dict[str, List[int]] = {}

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getNodeCoords(self, nodeIDs: Union[int, List[int], str]) -> array:
      """Return coordinates for one or more nodes as a NumPy array.

      Parameters
      ----------
      nodeIDs : Union[int, List[int], str]
        Single node ID (int), list of node IDs (list), or node group name (str).
        When a string is provided, it references a node group defined in the model.

      Returns
      -------
      numpy.array
        Array containing nodal coordinates. For a single node, returns a 1D array
        of shape (rank,). For multiple nodes, returns a 2D array of shape (n, rank)
        where n is the number of nodes and rank is the spatial dimension.

      Examples
      --------
      >>> nodes.getNodeCoords(5)           # Single node
      array([1.0, 2.0, 0.0])
      
      >>> nodes.getNodeCoords([1, 2, 3])   # Multiple nodes
      array([[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [2.0, 0.0, 0.0]])
      
      >>> nodes.getNodeCoords('Left')      # Node group
      array([[0.0, 0.0, 0.0],
             [0.0, 1.0, 0.0]])
      """

      if type(nodeIDs) == str:
        try:
            nodeIDs = self.groups[nodeIDs]
        except KeyError:
            logger.error(f"Node group '{nodeIDs}' not found.")
            sys.exit(1)

      return array(self.get(nodeIDs))
    
#
# 
#     

    def getNodeIDs(self, groupName: str) -> List[int]:
        """Return list of node IDs in the specified group.
    
        Parameters
        ----------
        groupName : str
            Name of the node group.
    
        Returns
        -------
        List[int]
            List of node IDs in the specified group.
    
        Examples
        --------
        >>> nodes.getNodeIDs('Left')
        [0, 1, 2]
        """
    
        try:
            return self.groups[groupName]
        except KeyError:
            logger.error(f"Node group '{groupName}' not found.")
            sys.exit(1)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def readFromFile(self, fname: str) -> None:
        """Read nodes and groups from legacy `.pro` file format.

        Parameters
        ----------
        fname : str
            Path to input file.
        """

        logger.info("  Reading nodes")
        logger.info("  -----------------------------------------------------------")

        fin = open(fname, 'r')

        line = fin.readline()

        while line:
            if line.replace(" ", "").startswith('<Nodes>'):
                self.readNodalCoords(fin)

            if line.replace(" ", "").startswith('gmsh'):
                ln = line.replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '').replace(';', '')
                ln = ln.split('=', 1)
                self.readGmshFile(ln[1][1:-1])
                break

            line = fin.readline()

        fin.close()

        fin = open(fname, 'r')

        line = fin.readline()

        while line:
            if line.replace(" ", "").startswith('<NodeGroup'):
                if 'name' in line:
                    label = (
                        line.split('=')[1]
                        .replace('\n', '')
                        .replace('>', '')
                        .replace(' ', '')
                        .replace('"', '')
                        .replace("'", '')
                    )
                    self.readNodegroup(fin, label)

            line = fin.readline()

        for key in self.groups:
            self.groups[key] = list(set(self.groups[key]))

        fin.close()
              
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def readGmshFile(self, fname: str) -> None:
        """Read nodes and groups from a Gmsh mesh file.

        Parameters
        ----------
        fname : str
            Path to the Gmsh file.
        """

        import meshio

        mesh = meshio.read(fname, file_format="gmsh")

        obj3d = ["pris", "pyra", "hexa", "wedg", "tetr"]

        self.rank = 2

        for key in mesh.cell_sets_dict:
            for typ in mesh.cell_sets_dict[key]:
                if (typ[:4] in obj3d):
                    self.rank = 3

        for nodeID, p in enumerate(mesh.points):
            self.add(nodeID, p[:self.rank])

        for key in mesh.cell_sets_dict:
            if key == "gmsh:bounding_entities":
                pass
            else:
                for typ in mesh.cell_sets_dict[key]:
                    for idx in mesh.cell_sets_dict[key][typ]:
                        iNodes = mesh.cells_dict[typ][idx]
                        for nodeID in iNodes:
                            self.addToGroup(key, nodeID)
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def addToGroup(self, modelType: str, ID: Union[int, str]) -> None:
      """Register a node id to a group.

      Parameters
      ----------
      modelType : str
        Group name.
      ID : Union[int, str]
        Node identifier (string parsed to int or already int).
      """

      if modelType not in self.groups:
        self.groups[modelType] = [int(ID)]
      else:
        self.groups[modelType].append(int(ID))
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __repr__(self) -> str:
        """Human-readable summary of node and group counts."""
        msg = f"  Number of nodes ............ {len(self):6d}\n"

        if len(self.groups) > 0:
            msg += f"    Number of  groups .......... {len(self.groups):6d}\n"
            msg += "    -----------------------------------\n"
            msg += "      name                       #nodes\n"
            msg += "      ---------------------------------\n"

            for name in self.groups:
                msg += f"      {name:<16s}           {len(self.groups[name]):6d} \n"

        return msg
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def readNodalCoords(self, fin: TextIO) -> None:
        """Read nodal coordinates from an open file stream.

        Parameters
        ----------
        fin : TextIO
            Open file handle positioned at the start of the `<Nodes>` block.
        """

        while True:
            line = fin.readline()

            if line.replace(" ", "").startswith('</Nodes>'):
                return

            line = re.sub(r'\s{2,}', ' ', line)
            a = line.split(';')

            for a in a[:-1]:
                b = a.strip().split(' ')

                if b[0].startswith("//") or b[0].startswith("#"):
                    break
                if len(b) > 1 and type(eval(b[0])) == int:
                    if self.rank == -1:
                        self.rank = len(b) - 1

                    self.add(eval(b[0]), [eval(crd) for crd in b[1:]])
          
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def readNodegroup(self, fin: TextIO, key: str) -> None:
        """Read a node group from an open file stream.

        Parameters
        ----------
        fin : TextIO
            Open file handle positioned at the start of a `<NodeGroup>` block.
        key : str
            Group name label.
        """

        while True:
            line = fin.readline()

            if line.replace(" ", "").startswith('</NodeGro'):
                return

            a = line.split()

            for b in a:
                if getType(b) == int:
                    self.addToGroup(key, b)

    def getRank(self) -> int:
        """Return the spatial dimension (rank) of the nodes.

        Returns
        -------
        int
            Spatial dimension (rank) of the nodes.
        """
        if self.rank == -1:
            self.rank = len(self.get(0))
            
        return self.rank
