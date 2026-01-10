# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import time
from typing import Any, Dict, Iterator, List, Optional, Sequence

import numpy as np
from numpy import zeros
from pyfem.util.logger import getLogger, separator

"""
Utilities for lightweight data structures used across PyFEM.

This module provides simple containers and helpers used by the solver and
pre/post-processing layers: cleaning input values, status tracking,
property containers, global data for assembly, and per-element storage.

All modifications in this file are non-functional documentation and type
annotations to improve readability and editor support; runtime behavior is
preserved.
"""

logger = getLogger()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def cleanVariable(a: Any) -> Any:
    """
    Convert a textual configuration value into a Python object.

    - 'true' -> True
    - 'false' -> False
    - otherwise: try `eval`, fall back to original string

    Args:
        a: Input value (typically a string) read from configuration.

    Returns:
        Parsed Python object or the original value if parsing fails.
    """
    if a == 'true':
        return True
    elif a == 'false':
        return False
    else:
        try:
            return eval(a)
        except Exception:
            return a

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class solverStatus:
    """Container for tracking solver progress and timing.

    Attributes:
        cycle (int): Current step/cycle number.
        iiter (int): Current iteration count within step.
        time (float): Current simulation time.
        time0 (float): Reference time (unused internally here).
        dtime (float): Time increment for a cycle.
        lam (float): Load or continuation parameter.
    """

    def __init__(self) -> None:
        self.cycle: int = 0
        self.iiter: int = 0
        self.time: float = 0.0
        self.time0: float = 0.0
        self.dtime: float = 0.0
        self.lam: float = 1.0

    def increaseStep(self) -> None:
        """Advance to the next solver step and reset iteration counter."""
        self.cycle += 1
        self.time += self.dtime
        self.iiter = 0
 
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
class Properties:
    """Simple attribute container built from a dictionary.

    Instances expose dictionary keys as attributes and support iteration over
    (name, value) pairs. This is a lightweight alternative to `types.SimpleNamespace`.
    """

    def __init__(self, dictionary: Dict[str, Any] = {}) -> None:
        for key in dictionary.keys():
            setattr(self, key, dictionary[key])

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __str__(self) -> str:
        """Return a multi-line representation listing public attributes."""
        myStr = ''
        for att in dir(self):
            # Ignore private members and standard routines
            if att.startswith('__'):
                continue
            myStr += 'Attribute: ' + att + '\n'
            myStr += str(getattr(self, att)) + '\n'
        return myStr

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __iter__(self) -> Iterator:
        """Iterate over (name, value) pairs for public attributes."""
        propsList: List = []
        for att in dir(self):
            # Ignore private members and standard routines
            if att.startswith('__'):
                continue
            propsList.append((att, getattr(self, att)))
        return iter(propsList)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def store(self, key: str, val: Any) -> None:
        """Store a property given a dotted key path or a simple key.

        Examples:
            store('alpha', 1)
            store('parent.child.value', 'x')
        """
        if '.' not in key:
            setattr(self, key, val)
        else:
            kets = key.split('.')
            props = self
            for y in kets[:-1]:
                props = getattr(props, y)
            setattr(props, kets[-1], cleanVariable(val))

#-------------------------------------------------------------------------------
# GlobalData - Global data container for finite element analysis
#-------------------------------------------------------------------------------   

class GlobalData(Properties):
    """Global data container for finite element analysis.
    
    This class holds the global state vectors, forces, boundary conditions,
    and manages I/O operations for nodal data. It inherits from Properties
    to provide flexible attribute storage.
    
    Attributes:
        nodes: Node container with coordinates and connectivity.
        elements: Element container with element definitions.
        dofs: Degree of freedom manager.
        state: Current state vector (displacements, temperatures, etc.).
        Dstate: Incremental state vector.
        fint: Internal force vector.
        fhat: Applied external force vector.
        velo: Velocity vector (for dynamic analysis).
        acce: Acceleration vector (for dynamic analysis).
        solverStatus: Solver status tracker.
        outputNames: List of output variable names.
    """
    
    def __init__(self, nodes: Any, elements: Any, dofs: Any) -> None:
        """Initialize the global data container.
        
        Args:
            nodes: Node container with nodal coordinates.
            elements: Element container with element definitions.
            dofs: Degree of freedom manager.
        """
        # Initialize base Properties class with nodes, elements, and dofs
        Properties.__init__(self, {'nodes': nodes, 'elements': elements, 'dofs': dofs})

        # Initialize state vectors
        self.state = zeros(len(self.dofs))
        self.Dstate = zeros(len(self.dofs))
        self.fint = zeros(len(self.dofs))
        self.fhat = zeros(len(self.dofs))

        # Initialize dynamic analysis vectors
        self.velo = zeros(len(self.dofs))
        self.acce = zeros(len(self.dofs))
        
        # Set solver status from elements
        self.solverStatus = elements.solverStat
       
        # Initialize output names list
        self.outputNames = []

    #---------------------------------------------------------------------------
    # readFromFile - Read external forces from input file
    #---------------------------------------------------------------------------

    def readFromFile(self, fname: str) -> None:
        """Read external forces from an input file.
        
        Parses the <ExternalForces> section of the input file and populates
        the fhat vector with prescribed external forces.
        
        Args:
            fname: Path to the input file containing external forces.
            
        Returns:
            None
            
        Note:
            Expected format in file:
            <ExternalForces>
            dofType[nodeID] = value;
            </ExternalForces>
        """
        # Open file and read all lines to check for ExternalForces section
        fin = open(fname)
        lines = fin.readlines()

        if not any('<ExternalForces>' in line for line in lines):
            logger.warning("No <ExternalForces> section found")
            return
        
        logger.info("Reading external forces ......")

        # Re-open file for parsing
        fin = open(fname)
      
        while True:
            line = fin.readline()
          
            # Find the start of ExternalForces section
            if line.startswith('<ExternalForces>') == True:
                while True:
                    line = fin.readline()

                    # Check for end of section
                    if line.startswith('</ExternalForces>') == True:
                        return
                
                    # Parse force specification line
                    a = line.strip().split(';')
                    
                    if len(a) == 2:
                        b = a[0].split('=')
                    
                        if len(b) == 2:
                            c = b[0].split('[')
                            
                            dofType = c[0]
                            nodeID = eval(c[1].split(']')[0])
                            
                            # Set external force value
                            self.fhat[self.dofs.getForType(nodeID, dofType)] = eval(b[1])
                  
        fin.close()

    #---------------------------------------------------------------------------
    # printNodes - Print nodal data to file or screen
    #---------------------------------------------------------------------------

    def printNodes(self, fileName: Optional[str] = None, inodes: Optional[List[int]] = None) -> None:
        """Print nodal data to file or screen.
        
        Outputs a formatted table containing nodal state values, internal forces,
        and any additional output variables.
        
        Args:
            fileName: Path to output file. If None, prints to screen.
            inodes: List of node IDs to print. If None, prints all nodes.
            
        Returns:
            None
        """
        # Determine output destination
        if fileName is None:
            f = None
        else:
            f = open(fileName, "w")

        # Get node list if not provided
        if inodes is None:
            inodes = list(self.nodes.keys())
        
        # Print header row
        print('     Node | ', file=f, end=' ')
        
        # Print DOF type headers
        for dofType in self.dofs.dofTypes:
            print(f"  {dofType:<10}", file=f, end=' ')

        # Print internal force headers if available
        if hasattr(self, 'fint'):
            for dofType in self.dofs.dofTypes:
                print(f" fint-{dofType:<6}", file=f, end=' ')

        # Print additional output variable headers
        for name in self.outputNames:
            print(f" {name:<11}", file=f, end=' ')

        print(" ", file=f)
        print('   ', ('-' * 81), file=f)

        # Print data for each node
        for nodeID in inodes:
            print(f'    {nodeID:4d}  | ', file=f, end=' ')
            
            # Print state values
            for dofType in self.dofs.dofTypes:
                print(f' {self.state[self.dofs.getForType(nodeID, dofType)]:10.3e} ', file=f, end=' ')
            
            # Print internal forces
            for dofType in self.dofs.dofTypes:
                print(f' {self.fint[self.dofs.getForType(nodeID, dofType)]:10.3e} ', file=f, end=' ')
            
            # Print additional output data
            for name in self.outputNames:
                print(f' {self.getData(name, nodeID):10.3e} ', file=f, end=' ')

            print(" ", file=f)
        print(" ", file=f)

        # Close file if writing to file
        if fileName is not None:
            f.close()

    #---------------------------------------------------------------------------
    # getData - Retrieve weighted output data for nodes
    #---------------------------------------------------------------------------
      
    def getData(self, outputName: str, inodes: Any) -> Any:
        """Retrieve weighted output data for specified nodes.
        
        Computes weighted average of output data for nodes. The weights are
        typically used to average data from multiple elements sharing a node.
        
        Args:
            outputName: Name of the output variable to retrieve.
            inodes: Node ID (int) or list of node IDs to retrieve data for.
            
        Returns:
            Weighted output value(s). Returns float for single node,
            list of values for multiple nodes.
        """
        # Get data and weight arrays
        data = getattr(self, outputName)
        weights = getattr(self, outputName + 'Weights')

        # Handle single node
        if type(inodes) is int:
            i = list(self.nodes.keys()).index(inodes)
            return data[i] / weights[i]
        else:
            # Handle multiple nodes
            outdata = []

            for row, w in zip(data[inodes], weights[inodes]):
                if w != 0:
                    outdata.append(row / w)
                else:
                    outdata.append(row)

            return outdata

    #---------------------------------------------------------------------------
    # resetNodalOutput - Clear all nodal output data
    #---------------------------------------------------------------------------

    def resetNodalOutput(self) -> None:
        """Clear all nodal output data and weights.
        
        Removes all output variables and their corresponding weights from
        the global data structure. Called at the start of each analysis step.
        
        Returns:
            None
        """
        # Remove all output data and weights
        for outputName in self.outputNames:
            delattr(self, outputName)
            delattr(self, outputName + 'Weights')

        # Clear output names list
        self.outputNames = []
    
    #---------------------------------------------------------------------------
    # close - Finalize analysis and print summary
    #---------------------------------------------------------------------------

    def close(self) -> None:
        """Finalize analysis and print execution summary.
        
        Prints the total elapsed time and a success message to the log.
        Called at the end of a successful analysis.
        
        Returns:
            None
        """
        from pyfem.util.plotUtils import plotTime
        
        logger.info("")
        separator("=")
        logger.info("  Total elapsed time.......... : " + 
                    plotTime(time.time() - self.startTime))
        logger.info("  PyFem analysis terminated successfully.")
        separator("=")
            
#-------------------------------------------------------------------------------
# elementData - Container for element-level data
#-------------------------------------------------------------------------------

class elementData():
    """Container for element-level state and computed quantities.
    
    This class stores element state vectors, element matrices (stiffness, mass),
    element force vectors, and output labels for post-processing.
    
    Attributes:
        state: Element state vector (displacements, etc.).
        Dstate: Element incremental state vector.
        stiff: Element tangent stiffness matrix.
        fint: Element internal force vector.
        mass: Element mass matrix.
        lumped: Element lumped mass vector.
        diss: Element dissipation (energy dissipated).
        outlabel: List of output variable labels for this element.
    """

    def __init__(self, elstate: np.ndarray, elDstate: np.ndarray) -> None:
        """Initialize element data container.
        
        Args:
            elstate: Element state vector.
            elDstate: Element incremental state vector.
        """
        # Get number of DOFs from state vector
        nDof = len(elstate)

        # Store state vectors
        self.state = elstate
        self.Dstate = elDstate
        
        # Initialize element matrices and vectors
        self.stiff = zeros(shape=(nDof, nDof))
        self.fint = zeros(shape=(nDof))
        self.mass = zeros(shape=(nDof, nDof))
        self.lumped = zeros(shape=(nDof))
        self.diss = 0.0

        # Initialize output labels
        self.outlabel = []
   
    def __str__(self) -> str:
        """Return string representation of element state.
        
        Returns:
            String representation of the state vector.
        """
        return str(self.state)
     
