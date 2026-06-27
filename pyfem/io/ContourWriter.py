# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Any, List
from pyfem.util.BaseModule import BaseModule


class ContourWriter(BaseModule):
    """Output module for writing contour data along specified node sets.
    
    Writes nodal coordinates, degrees of freedom, and output variables along
    a contour line at specified intervals during the simulation.
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the ContourWriter.
        
        Args:
            props: Properties dictionary containing module configuration.
            globdat: Global data object containing simulation state.
        """
        self.prefix = globdat.prefix
        self.interval = 1

        BaseModule.__init__(self, props)
        
        self.k = 0
        self.columndata: List[Any] = []

    def run(self, props: Any, globdat: Any) -> None:
        """Write contour data for the current simulation step.
        
        Outputs nodal coordinates, DOF values, and element output variables
        for all nodes in the specified contour at the current cycle.
        
        Args:
            props: Properties dictionary (not used in this method).
            globdat: Global data object containing simulation state.
        """
        if not globdat.solverStatus.cycle % self.interval == 0:
            return
     
        self.writeHeader()
 
        crd = globdat.nodes.getNodeCoords(self.nodes[0])
        outfile = open(self.prefix + '-contour-' + str(self.k) + '.out', 'w')
        
        outfile.write(f'#Node  {"x-coor":<10} {"y-coor":<10}')
  
        if len(crd) == 3:
            outfile.write(f'{"z-coor":<10} ')

        for dofType in globdat.dofs.dofTypes:
            outfile.write(f'{dofType:<10} ')

        for name in globdat.outputNames:
            outfile.write(f'{name:<10} ')

        outfile.write('\n')

        for iNod in self.nodes:
            crd = globdat.nodes.getNodeCoords(iNod)
            outfile.write(f'{iNod:4d} {crd[0]:10.3e} {crd[1]:10.3e}')

            if len(crd) == 3:
                outfile.write(f' {crd[2]:10.3e}')

            for dofType in globdat.dofs.dofTypes:
                outfile.write(f' {globdat.state[globdat.dofs.getForType(iNod, dofType)]:10.3e}')

            for name in globdat.outputNames:
                stress = globdat.getData(name, list(range(len(globdat.nodes))))
                outfile.write(f' {stress[iNod]:10.3e}')

            outfile.write('\n')
                
        outfile.close()
  
        self.k = self.k + 1
