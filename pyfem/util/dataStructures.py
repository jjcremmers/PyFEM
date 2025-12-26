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

import time
from typing import Any, Dict, Iterator, List, Optional, Sequence

import numpy as np
from numpy import zeros
from pyfem.util.logger import getLogger

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
#
#-------------------------------------------------------------------------------   

class GlobalData ( Properties ):
  
  def __init__ ( self, nodes, elements, dofs ):

    Properties.__init__( self, { 'nodes' : nodes, 'elements' : elements, 'dofs' : dofs } )

    self.state  = zeros( len( self.dofs ) )
    self.Dstate = zeros( len( self.dofs ) )
    self.fint   = zeros( len( self.dofs ) )
    self.fhat   = zeros( len( self.dofs ) )

    self.velo   = zeros( len( self.dofs ) )
    self.acce   = zeros( len( self.dofs ) )
    
    self.solverStatus = elements.solverStat
   
    self.outputNames = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def readFromFile( self , fname ):

    logger.info("Reading external forces ......")

    fin = open( fname )
  
    while True:
      line = fin.readline()  
  
      if line.startswith('<ExternalForces>') == True:
        while True:
          line = fin.readline()  

          if line.startswith('</ExternalForces>') == True:
            return
        
          a = line.strip().split(';')
          
          if len(a) == 2:
            b = a[0].split('=')
        
            if len(b) == 2:
              c = b[0].split('[')
              
              dofType = c[0]
              nodeID  = eval(c[1].split(']')[0])
              
              self.fhat[self.dofs.getForType(nodeID,dofType)] = eval(b[1])
              
    fin.close()

#---------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------

  def printNodes( self , fileName = None , inodes=None ):

    if fileName is None:
      f = None
    else:
      f = open(fileName,"w")

    if inodes is None:
      inodes = list(self.nodes.keys()) 
	
    print('   Node | ', file = f , end=' ')
    
    for dofType in self.dofs.dofTypes:
      print(f"  {dofType:<10}", file=f, end=' ')

    if hasattr( self , 'fint' ):
      for dofType in self.dofs.dofTypes:
        print(f" fint-{dofType:<6}", file=f, end=' ')

    for name in self.outputNames:
      print(f" {name:<11}", file=f, end=' ')

    print(" ", file = f ) 
    print(('-' * 100), file = f )

    for nodeID in inodes:
      print(f'  {nodeID:4d}  | ', file=f, end=' ')
      for dofType in self.dofs.dofTypes:
        print(f' {self.state[self.dofs.getForType(nodeID, dofType)]:10.3e} ', file=f, end=' ')
      for dofType in self.dofs.dofTypes:
        print(f' {self.fint[self.dofs.getForType(nodeID, dofType)]:10.3e} ', file=f, end=' ')
      for name in self.outputNames:
        print(f' {self.getData(name, nodeID):10.3e} ', file=f, end=' ')

      print(" ", file = f )
    print(" ", file = f )

    if fileName is not None:
      f.close()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
      
  def getData( self , outputName , inodes ):

    data    = getattr( self, outputName )
    weights = getattr( self, outputName + 'Weights' )

    if type(inodes) is int:
      i = list(self.nodes.keys()).index(inodes)
      return data[i] / weights[i]
    else:
      outdata=[]

      for row,w in zip(data[inodes],weights[inodes]):
        if w != 0:
          outdata.append(row/w)
        else:
          outdata.append(row)

      return outdata

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def resetNodalOutput ( self ):

    for outputName in self.outputNames:
      delattr( self , outputName )
      delattr( self , outputName + 'Weights' )

    self.outputNames=[]
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def close ( self ):
 
    from pyfem.util.plotUtils   import plotTime 
    
    logger.info("")
    logger.info("=============================================================")
    logger.info("  Total elapsed time.......... : " + 
                   plotTime(time.time()-self.startTime))
    logger.info("  PyFem analysis terminated successfully.")
    logger.info("=============================================================")  
            
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class elementData():

  def __init__( self , elstate , elDstate ):

    nDof        = len(elstate)

    self.state  = elstate
    self.Dstate = elDstate
    self.stiff  = zeros( shape=( nDof,nDof ) )
    self.fint   = zeros( shape=(nDof) )
    self.mass   = zeros( shape=( nDof,nDof ) )
    self.lumped = zeros( shape=(nDof) )
    self.diss   = 0.0

    self.outlabel = []
   
  def __str__( self ):

    return self.state
     
