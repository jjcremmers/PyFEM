# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

from typing import Any, List

import h5py
import numpy as np

from pyfem.util.BaseModule import BaseModule


class ROMSnapshotWriter(BaseModule):
    """Collect full-order snapshots for reduced-order modeling workflows.

    The writer appends state vectors to an HDF5 file. On the first call it
    also stores mesh topology, nodal coordinates, and a node-wise displacement
    DOF map so the ROM preparation and online stages can reconstruct the model
    layout from the same database.
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the ROM snapshot writer.

        Args:
            props: Module configuration properties.
            globdat: Global data object providing the output prefix and mesh
                metadata.
        """
        self.prefix = globdat.prefix
        self.extension = ".h5"
        self.dispDofs = ["u", "v", "w"]

        BaseModule.__init__(self, props)

        if not hasattr(self, "filename"):
            self.filename = self.prefix + self.extension

        self.h5file = h5py.File(self.filename, "a")

    def run(self, props: Any, globdat: Any) -> None:
        """Store the current state vector in the ROM HDF5 file.

        Args:
            props: Module configuration properties. Not used directly.
            globdat: Global data object containing the current state vector,
                nodes, elements, and DOF numbering.
        """
        state_size = len(globdat.state)

        if "state" not in self.h5file:
            self.h5file.create_dataset(
                "state",
                maxshape=(None, state_size),
                chunks=True,
                data=globdat.state.reshape(1, state_size),
                dtype="f",
            )
            self._write_mesh_data(globdat)
        else:
            state_data = self.h5file["state"]
            state_data.resize((state_data.shape[0] + 1), axis=0)
            state_data[-1, :] = globdat.state

    def _write_mesh_data(self, globdat: Any) -> None:
        """Write mesh and nodal metadata to the HDF5 file.

        Args:
            globdat: Global data object containing mesh topology and nodal
                information.
        """
        self.h5file.create_group("elements")
        self.h5file.create_group("nodes")

        self._write_element_data(globdat)
        self._write_node_data(globdat)

    def _write_element_data(self, globdat: Any) -> None:
        """Write element offsets and connectivity.

        Args:
            globdat: Global data object containing the element set.
        """
        offsets: List[int] = []
        connectivity: List[int] = []

        offset = 0
        for element in globdat.elements:
            offset += len(element)
            offsets.append(offset)
            connectivity.extend(element)

        connectivity_array = np.array(
            globdat.nodes.getIndices(connectivity), dtype=int
        )
        offsets_array = np.array(offsets, dtype=int)

        self.h5file["elements"].create_dataset(
            "offsets", offsets_array.shape, dtype="i", data=offsets_array
        )
        self.h5file["elements"].create_dataset(
            "connectivity",
            connectivity_array.shape,
            dtype="i",
            data=connectivity_array,
        )

    def _write_node_data(self, globdat: Any) -> None:
        """Write nodal coordinates and displacement-DOF indices.

        Args:
            globdat: Global data object containing nodes and DOF numbering.
        """
        coordinates = []

        for node_id in list(globdat.nodes.keys()):
            coordinates.append(globdat.nodes.getNodeCoords(node_id))

        coordinates_array = np.array(coordinates, dtype=float)
        dof_map = self._build_dof_map(globdat, coordinates_array.shape[1])

        self.h5file["nodes"].create_dataset(
            "coordinates",
            coordinates_array.shape,
            dtype="f",
            data=coordinates_array,
        )
        self.h5file.create_dataset("dofs", dof_map.shape, dtype="i", data=dof_map)

    def _build_dof_map(self, globdat: Any, rank: int) -> np.ndarray:
        """Build a node-wise map from displacement DOFs to global indices.

        Args:
            globdat: Global data object containing the degree-of-freedom space.
            rank: Spatial dimension inferred from the nodal coordinates.

        Returns:
            ndarray: Integer array of shape ``(nnode, rank)`` with the global
                displacement DOF indices for each node.
        """
        active_disp_dofs = self.dispDofs[:rank]
        dof_map = []

        for node_id in list(globdat.nodes.keys()):
            node_dofs = []
            for disp_dof in active_disp_dofs:
                if disp_dof in globdat.dofs.dofTypes:
                    node_dofs.append(globdat.dofs.getForType(node_id, disp_dof))
                else:
                    node_dofs.append(-1)
            dof_map.append(node_dofs)

        return np.array(dof_map, dtype=int)