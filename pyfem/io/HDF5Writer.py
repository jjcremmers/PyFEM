# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Any, List
import h5py
import numpy as np

from pyfem.util.BaseModule import BaseModule


class HDF5Writer(BaseModule):
    """Module for writing simulation data to HDF5 format.
    
    Exports node coordinates, element connectivity, displacements, and
    other output fields to HDF5 files for post-processing and visualization.
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the HDF5Writer module.
        
        Args:
            props: Properties dictionary containing module configuration.
            globdat: Global data object containing simulation state.
        """
        self.prefix = globdat.prefix
        self.extension = ".h5"
        
        self.dispDofs = ["u", "v", "w"]
        self.extraFields = ["rx", "ry", "rz", "temp", "pres", "phase"]
        self.singleFile = True
        
        BaseModule.__init__(self, props)
        
        if not hasattr(props, "interval"):
            self.interval = 1
        
        if self.singleFile:
            f = h5py.File(self.prefix + self.extension, "w")
            f.attrs["cycleCount"] = 0

    def run(self, props: Any, globdat: Any) -> None:
        """Write simulation data to HDF5 file.
        
        Writes mesh data, displacements, and output variables at specified
        intervals. Handles both single-file and multi-file modes, as well
        as eigenmode output.
        
        Args:
            props: Properties dictionary (not used in this method).
            globdat: Global data object containing current simulation state.
        """
        cycle = globdat.solverStatus.cycle
        
        if hasattr(globdat, "eigenvecs"):
            
            self.writeHeader()
            name = str(self.prefix + self.extension)
            
            f = h5py.File(name, "w")
            
            self.writeCycle(f, globdat, method="modes")
        
        elif cycle % self.interval == 0:
            
            self.writeHeader()
            
            if self.singleFile:
                f = h5py.File(self.prefix + self.extension, "a")
                
                self.cycle = f.attrs["cycleCount"]
                
                self.cycle += 1
                
                f.attrs["cycleCount"] = self.cycle
                
                gName = "cycle" + str(self.cycle)
                f.create_group(gName)
                
                self.writeCycle(f[gName], globdat)
                
            else:
                name = str(self.prefix + "_" + str(self.cycle) + self.extension)
            
                f = h5py.File(name, "w")
                    
                self.writeCycle(f, globdat)
                
                self.cycle += 1

    def writeCycle(self, cdat: Any, globdat: Any, method: str = "all") -> None:
        """Write data for a single simulation cycle to HDF5.
        
        Writes mesh topology, node coordinates, displacements, and output
        variables for the current cycle or eigenmode data.
        
        Args:
            cdat: HDF5 group or file object to write data to.
            globdat: Global data object containing simulation state.
            method: Output method - \"all\" for full data or \"modes\" for eigenmodes.
        """
        cdat.create_group("elements")
        
        elemCount = []
        connectivity = []
        
        i0 = 0
        for elem in globdat.elements:
            i0 = i0 + len(elem)
            elemCount.append(i0)
            connectivity.extend(elem)
        
        connectivity = np.array(globdat.nodes.getIndices(connectivity), dtype=int)
        elemCount = np.array(elemCount, dtype=int)
        elemIDs = np.array(globdat.elements.getIndices(), dtype=int)
        familyIDs = np.array(globdat.elements.getFamilyIDs(), dtype=int)
        
        cdat["elements"].create_dataset("offsets", elemCount.shape,
                                        dtype='i', data=elemCount)
        cdat["elements"].create_dataset("connectivity", connectivity.shape,
                                        dtype='i', data=connectivity)
        cdat["elements"].create_dataset("elementIDs", elemIDs.shape,
                                        dtype='i', data=elemIDs)
        cdat["elements"].create_dataset("familyIDs", familyIDs.shape,
                                        dtype='i', data=familyIDs)
        
        cdat.create_group("elementGroups")
        
        for key in globdat.elements.groups:
            elementIDs = np.array(globdat.elements.getIndices(globdat.elements.groups[key]), dtype=int)
            cdat["elementGroups"].create_dataset(key, elementIDs.shape, dtype='i', data=elementIDs)
        
        cdat.create_group("nodes")
        
        coordinates = []
        
        for nodeID in list(globdat.nodes.keys()):
            coordinates.append(globdat.nodes.getNodeCoords(nodeID))
        
        coordinates = np.array(coordinates, dtype=float)
        nodeIDs = np.array(globdat.nodes.getIndices(), dtype=int)
        
        cdat["nodes"].create_dataset("coordinates", coordinates.shape,
                                     dtype='f', data=coordinates)
        cdat["nodes"].create_dataset("nodeIDs", nodeIDs.shape,
                                     dtype='i', data=nodeIDs)
        
        dofs = self.dispDofs[:coordinates.shape[1]]
        
        cdat.create_group("nodeGroups")
        
        for key in globdat.nodes.groups:
            nodeIDs = np.array(globdat.nodes.getIndices(globdat.nodes.groups[key]), dtype=int)
            cdat["nodeGroups"].create_dataset(key, nodeIDs.shape, dtype='i', data=nodeIDs)
        
        cdat.create_group("nodeData")
        
        if method == "all":
            displacements = []

            for nodeID in list(globdat.nodes.keys()):
                d = []
                for dispDof in dofs:
                    if dispDof in globdat.dofs.dofTypes:
                        d.append(globdat.state[globdat.dofs.getForType(nodeID, dispDof)])
                    else:
                        d.append(0.)
                displacements.append(d)
            
            displacements = np.array(displacements, dtype=float)
    
            cdat["nodeData"].create_dataset("displacements", displacements.shape, dtype='f', data=displacements)
            
            for field in self.extraFields:
                if field in globdat.dofs.dofTypes:
                    output = []
                        
                    for nodeID in list(globdat.nodes.keys()):
                        output.append(globdat.state[globdat.dofs.getForType(nodeID, field)])
        
                    output = np.array(output, dtype=float)
            
                    cdat["nodeData"].create_dataset(field, output.shape, dtype='f', data=output)
            
            for name in globdat.outputNames:
                output = globdat.getData(name, list(range(len(globdat.nodes))))
                
                output = np.array(output, dtype=float)
            
                cdat["nodeData"].create_dataset(name, output.shape, dtype='f', data=output)
    
            if hasattr(globdat, "elementData"):
                cdat.create_group("elemData")
                elemData = globdat.elementData
                            
                for name in elemData.outputNames:
                    data = getattr(elemData, name)
                
                    output = np.array(data, dtype=float)
            
                    cdat["elemData"].create_dataset(name, output.shape, dtype='f', data=output)
        
        elif method == "modes":
        
            allModes = []
            
            for iMod, eigenvecs in enumerate(globdat.eigenvecs.T):
                mode = []
            
                for nodeID in list(globdat.nodes.keys()):
                    d = []
                    for dispDof in dofs:
                        if dispDof in globdat.dofs.dofTypes:
                            d.append(globdat.eigenvecs[globdat.dofs.getForType(nodeID, dispDof)])
                        else:
                            d.append(0.)
                    mode.append(d)
                
                allModes.append(mode)
            
            allModes = np.array(allModes, dtype=float)
    
            cdat["nodeData"].create_dataset("modes", allModes.shape, dtype='f', data=allModes)
            
            eigenvals = globdat.eigenvals
            
            cdat.create_dataset("eigenvals", eigenvals.shape, dtype='f', data=eigenvals)
