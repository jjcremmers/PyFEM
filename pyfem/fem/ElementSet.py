# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import re
from importlib import import_module
from typing import Any, Dict, Iterator, List, Union
from pyfem.util.itemList import itemList
from pyfem.util.logger import getLogger, logHeader, separator, logVariable
from pyfem.util.dataStructures import solverStatus, Properties

logger = getLogger()


class ElementSet(itemList):
    """Container for finite elements and element groups.
    
    Manages a collection of finite element objects indexed by element ID.
    Elements are organized into named groups based on their model type (material,
    element type, etc.). Provides methods for iterating over elements, reading
    from input files, and managing element groups.
    
    Attributes:
        nodes: NodeSet object containing the mesh nodes.
        props: Properties object containing model and material properties.
        solverStat: Solver status tracker for time stepping and iterations.
        groups: Dictionary mapping group names to lists of element IDs.
    
    Examples:
        >>> elements = ElementSet(nodes, props)
        >>> elements.readFromFile('input.pro')
        >>> for element in elements:
        ...     # Process each element
    """

    def __init__(self, nodes: Any, props: Properties) -> None:
        """Initialize an empty ElementSet.
        
        Args:
            nodes: NodeSet object containing mesh nodes.
            props: Properties object with model and material definitions.
        """
        itemList.__init__(self)
        
        self.nodes = nodes
        self.props = props
        self.solverStat = solverStatus()
        self.groups: Dict[str, List[int]] = {}

    #---------------------------------------------------------------------------
    # __iter__ - Iterate over all elements
    #---------------------------------------------------------------------------

    def __iter__(self) -> Iterator:
        """Iterate over all elements in all groups.
        
        Returns:
            Iterator over all element objects in the set.
        """
        elements = []

        for groupName in self.iterGroupNames():
            for element in self.iterElementGroup(groupName):
                elements.append(element)
           
        return iter(elements)

    #---------------------------------------------------------------------------
    # __repr__ - Return string representation
    #---------------------------------------------------------------------------
        
    def __repr__(self) -> str:
        """Return human-readable summary of element and group counts.
        
        Returns:
            Formatted string showing total element count and group information.
        """
        msg = f"  Number of elements ......... {len(self):6d}\n"
        
        if len(self.groups) > 0:
            msg += f"    Number of  groups .......... {len(self.groups):6d}\n"
            msg += "    -----------------------------------\n"
            msg += "      name                       #elems\n"
            msg += "      ---------------------------------\n"
            for name in self.groups:
                msg += f"      {name:<16s}           {len(self.groups[name]):6d}\n"
        
        return msg

    #---------------------------------------------------------------------------
    #
    #---------------------------------------------------------------------------

    def logInfo(self) -> None:
        """Log element and group information using the logger.
        
        Outputs the same information as __repr__ but using the logger,
        with properly formatted output for number of elements and groups.
        """
        logVariable("Number of elements", len(self))
        
        if len(self.groups) > 0:
            logVariable("Number of groups", len(self.groups))
            logger.info("")
            logHeader("    name",  "#elems")

            for name in self.groups:
                logVariable(("    " + name),len(self.groups[name])) 

    #---------------------------------------------------------------------------
    # getDofTypes - Get all unique DOF types from elements
    #---------------------------------------------------------------------------

    def getDofTypes(self) -> List[str]:
        """Get all unique degree of freedom (DOF) types from all elements.
        
        Iterates through all elements and collects unique DOF types (e.g., 'u', 'v',
        'w' for displacements, 'temp' for temperature, etc.).
        
        Returns:
            List of unique DOF type strings.
        """
        dofTypes = []

        for element in self:
            for dofType in element.dofTypes:
                if dofType not in dofTypes:
                    dofTypes.append(dofType)

        return dofTypes
        
    #---------------------------------------------------------------------------
    # readFromFile - Read elements from input file
    #---------------------------------------------------------------------------
        
    def readFromFile(self, fname: str) -> None:
        """Read element definitions from an input file.
        
        Parses elements from a legacy .pro file format or references a Gmsh file.
        Supports <Elements> blocks with element definitions or gmsh file references.
        
        Args:
            fname: Path to input file (.pro format or referencing Gmsh file).
            
        Notes:
            Format for <Elements> block:
            elemID modelName nodeID1 nodeID2 ...; 
            
            Format for Gmsh reference:
            gmsh = "filename.msh";
        """
        logger.info("Reading elements")
        separator()

        fin = open(fname)
      
        while True:
            line = fin.readline()
          
            if line.startswith('<Elements>') == True:
                while True:
                    line = fin.readline()

                    if line.startswith('</Elements>') == True:
                        return
                
                    line = re.sub(r'\s{2,}', ' ', line)
                    a = line.split(';')
         
                    for a0 in a[:-1]:
                        b = a0.strip().split(' ')
                        
                        if b[0].startswith("//") or b[0].startswith("#"):
                            break
                        if len(b) > 1 and type(eval(b[0])) == int:
                            self.add(eval(b[0]), eval(b[1]), [eval(nodeID) for nodeID in b[2:]])

            elif line.startswith('gmsh') == True:
                ln = line.replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '').replace(';', '')
                ln = ln.split('=', 1)
                self.readGmshFile(ln[1][1:-1])
                return
                
        fin.close()
            
    #---------------------------------------------------------------------------
    # readGmshFile - Read elements from Gmsh mesh file
    #---------------------------------------------------------------------------

    def readGmshFile(self, fname: str) -> None:
        """Read elements from a Gmsh mesh file.
        
        Uses meshio to parse Gmsh format and extract elements. Physical groups
        from Gmsh become element groups in PyFEM.
        
        Args:
            fname: Path to the Gmsh file (.msh format).
            
        Requires:
            meshio package for reading Gmsh files.
        """
        import meshio
        
        mesh = meshio.read(fname, file_format="gmsh")
        
        elemID = 0

        for key in mesh.cell_sets_dict:
            for typ in mesh.cell_sets_dict[key]:
                for idx in mesh.cell_sets_dict[key][typ]:
                    iNodes = mesh.cells_dict[typ][idx]
                    self.add(elemID, key, iNodes.tolist())
                    elemID = elemID + 1
                                     
    #---------------------------------------------------------------------------
    # add - Add an element to the set
    #---------------------------------------------------------------------------

    def add(self, ID: int, modelName: str, elementNodes: List[int]) -> None:
        """Add a finite element to the element set.
        
        Creates an element instance of the specified type with given nodes,
        validates node IDs, and adds the element to the appropriate group.
        
        Args:
            ID: Unique element identifier.
            modelName: Name of the model/material (must exist in props).
            elementNodes: List of node IDs that define the element connectivity.
            
        Raises:
            RuntimeError: If model is missing 'type' attribute or node ID invalid.
            ImportError: If element module or class cannot be found.
        """
        # Check if the model exists
        
        if hasattr(self.props, modelName):
        
            modelProps = getattr(self.props, modelName)

            # Check if the model has a type
            if not hasattr(modelProps, 'type'):
                raise RuntimeError('Missing type for model ' + modelName)
            
            modelType = getattr(modelProps, 'type')
     
            modelProps.rank = self.nodes.rank
            modelProps.solverStat = self.solverStat

            try:
                model = import_module(f"pyfem.elements.{modelType}")
            except ModuleNotFoundError as e:
                raise ImportError(
                    f"Solver module 'pyfem.elements.{modelType}' not found. "
                    f"Check the 'type' in your input file."
                ) from e

            try:
                model_cls = getattr(model, modelType)
            except AttributeError as e:
                raise ImportError(
                    f"Class '{modelType}' not found in module 'pyfem.elements.{modelType}'. "
                    f"Ensure the class name matches the file name."
                ) from e
     
            elem = model_cls(elementNodes, modelProps)

            # Check if the node IDs are valid:

            for nodeID in elem.getNodes():
                if not nodeID in self.nodes:
                    raise RuntimeError('Node ID ' + str(nodeID) + ' does not exist')

            # Add the element to the element set:

            itemList.add(self, ID, elem)

            # Add the element to the correct group:

            self.addToGroup(modelName, ID)

    #---------------------------------------------------------------------------
    # addToGroup - Add element to a named group
    #---------------------------------------------------------------------------

    def addToGroup(self, modelType: str, ID: int) -> None:
        """Register an element ID to a named group.
        
        Creates the group if it doesn't exist, otherwise appends the element ID.
        
        Args:
            modelType: Name of the element group.
            ID: Element identifier to add to the group.
        """
        if modelType not in self.groups:
            self.groups[modelType] = [ID]
        else:
            self.groups[modelType].append(ID)

    #---------------------------------------------------------------------------
    # addGroup - Create a new group with specified elements
    #---------------------------------------------------------------------------

    def addGroup(self, groupName: str, groupIDs: List[int]) -> None:
        """Create or replace a named group with specified element IDs.
        
        Args:
            groupName: Name for the element group.
            groupIDs: List of element IDs to include in the group.
        """
        self.groups[groupName] = groupIDs

    #---------------------------------------------------------------------------
    # iterGroupNames - Get group names
    #---------------------------------------------------------------------------

    def iterGroupNames(self) -> Dict[str, List[int]]:
        """Return the dictionary of group names and their element IDs.
        
        Returns:
            Dictionary mapping group names to lists of element IDs.
        """
        return self.groups

    #---------------------------------------------------------------------------
    # iterElementGroup - Iterate over elements in a group
    #---------------------------------------------------------------------------

    def iterElementGroup(self, groupName: Union[str, List[str]]) -> Iterator:
        """Iterate over element objects in specified group(s).
        
        Args:
            groupName: Name of group to iterate, 'All' for all elements,
                      or list of group names to iterate multiple groups.
                      
        Returns:
            Iterator over element objects in the specified group(s).
            
        Examples:
            >>> for elem in elements.iterElementGroup('Material1'):
            ...     # Process elements in Material1 group
            >>> for elem in elements.iterElementGroup(['Mat1', 'Mat2']):
            ...     # Process elements in both groups
        """
        if groupName == "All":
            return iter(self)
        elif isinstance(groupName, list):
            elems = []
            for name in groupName:
                elems += self.get(self.groups[name])
            return iter(elems)
        else:
            return iter(self.get(self.groups[groupName]))

    #---------------------------------------------------------------------------
    # elementGroupCount - Get number of elements in group
    #---------------------------------------------------------------------------

    def elementGroupCount(self, groupName: Union[str, List[str]]) -> int:
        """Return the number of elements in specified group(s).
        
        Args:
            groupName: Name of group, 'All' for all elements,
                      or list of group names.
                      
        Returns:
            Total count of elements in the specified group(s).
        """

        if groupName == "All":
            return len(self)
        elif isinstance(groupName, list):
            length = 0
            for name in groupName:
                length += len(self.groups[name])
            return length
        else:
            return len(self.groups[groupName])
          
    #---------------------------------------------------------------------------
    # getFamilyIDs - Get element family indices
    #---------------------------------------------------------------------------

    def getFamilyIDs(self) -> List[int]:
        """Get family type indices for all elements.
        
        Returns a list of integer indices representing the element family
        (CONTINUUM=0, INTERFACE=1, SURFACE=2, BEAM=3, SHELL=4) for each element.
        
        Returns:
            List of family indices, one per element.
        """
        familyIDs = []
        fam = ["CONTINUUM", "INTERFACE", "SURFACE", "BEAM", "SHELL"]
        
        for elem in self:
            
            familyIDs.append(fam.index(elem.family))
            
        return familyIDs

    #---------------------------------------------------------------------------
    # commitHistory - Commit element history variables
    #---------------------------------------------------------------------------

    def commitHistory(self) -> None:
        """Commit history variables for all elements.
        
        Calls commitHistory() on all element objects to update internal
        state variables after a successful time/load step. Typically used
        for storing plastic strains, damage variables, etc.
        """
        for element in list(self.values()):
            element.commitHistory()
