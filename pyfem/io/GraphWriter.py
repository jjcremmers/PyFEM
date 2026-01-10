# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Any, List
import matplotlib.pyplot as plt
from numpy import ndarray, zeros

from pyfem.util.BaseModule import BaseModule
from pyfem.util.dataStructures import Properties


class GraphWriter(BaseModule):
    """Module for writing graph data to files and optionally displaying on screen.
    
    Collects specified output variables at each time step and writes them to
    a file. Can also display real-time plots using matplotlib.
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the GraphWriter module.
        
        Args:
            props: Properties dictionary containing module configuration,
                   including columns to output.
            globdat: Global data object containing simulation state.
        """
        self.prefix = globdat.prefix
        self.extension = ".out"
        self.onScreen = False

        BaseModule.__init__(self, props)

        if not hasattr(self, "filename"):
            self.filename = self.prefix + self.extension
        
        self.columndata: List[Any] = []

        for i, col in enumerate(self.columns):

            if hasattr(self, col):
                colProps = getattr(self, col)
            else:
                colProps = Properties()
                
            if not hasattr(colProps, "type"):
                colProps.type = col
            
            if not hasattr(colProps, "factor"):
                colProps.factor = 1.0
                
            if hasattr(colProps, "node"):
                if type(colProps.node) == str:
                    colProps.node = globdat.nodes.groups[colProps.node]

            self.columndata.append(colProps)

        if self.onScreen:
            globdat.onScreen = True

            self.fig = plt.figure(figsize=(3, 4), dpi=160)
            self.ax1 = plt.subplot()
            
        self.outfile = open(self.filename, 'w')

        if self.onScreen:
            self.output: List[List[float]] = []
            
        self.outfile = open(self.filename, 'w')

        self.run(props, globdat)

    def run(self, props: Any, globdat: Any) -> None:
        """Write graph data for the current simulation step.
        
        Collects data from all configured columns and writes them to the
        output file. If onScreen is enabled, updates the matplotlib plot.
        
        Args:
            props: Properties dictionary (not used in this method).
            globdat: Global data object containing current simulation state.
        """
        self.writeHeader()
        
        a: List[float] = []

        for i, col in enumerate(self.columndata):
            
            if col.type in globdat.outputNames:
                data = globdat.getData(col.type, col.node)
                
            elif hasattr(globdat, col.type):
                b = getattr(globdat, col.type)
                if type(b) is ndarray:
                    if type(col.node) is list:
                        data = 0.0
                        for nod in col.node:
                            data += b[globdat.dofs.getForType(int(nod), col.dof)]
                    else:
                        data = b[globdat.dofs.getForType(col.node, col.dof)]
                else:
                    data = b
                    
            elif col.type in globdat.outputNames:
                data = globdat.getData(col.type, col.node)
                
            elif hasattr(globdat.solverStatus, col.type):
                data = getattr(globdat.solverStatus, col.type)
                
            else:
                data = 0.0
        
            data = data * col.factor

            a.append(data)
        
            self.outfile.write(str(data) + ' ',)
            self.outfile.flush()

        self.outfile.write('\\n')

        if self.onScreen:
            self.output.append(a)
            
            plt.sca(self.ax1)
            plt.cla()
            
            plt.xlabel(self.columns[0])
            plt.ylabel(self.columns[1])
            
            plt.plot([x[0] for x in self.output], [x[1] for x in self.output], 'ro-')
        
            plt.pause(0.001)
            
            self.fig.savefig(self.prefix + '.png')
        
        if not globdat.active:
            self.outfile.close
