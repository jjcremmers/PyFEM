# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

from typing import Any

import h5py
import numpy as np
from scipy.linalg import svd

from pyfem.util.BaseModule import BaseModule
from pyfem.util.logger import getLogger

logger = getLogger()


class ROMBasisBuilder(BaseModule):
    """Prepare reduced-order bases from collected full-order snapshots.

    The solver reads an HDF5 snapshot database, applies the requested ROM basis
    construction method, and stores the resulting modes and singular values in
    the same file. At present only POD-based preparation is implemented.
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the ROM basis builder module.

        Args:
            props: Solver configuration properties.
            globdat: Global data object. It is accepted for interface
                consistency with the other solver modules.
        """
        self.method = "POD"

        BaseModule.__init__(self, props)

    def run(self, props: Any, globdat: Any) -> None:
        """Build reduced-order data from the configured snapshot database.

        Args:
            props: Solver configuration properties.
            globdat: Global data object whose ``active`` flag is cleared after
                basis construction finishes.

        Raises:
            NotImplementedError: If the requested preparation method is not
                available.
        """
        self.writeHeader()

        logger.info("Preparing reduced-order basis ........")
        logger.info("    snapshot file    : %s", self.filename)
        logger.info("    method           : %s", self.method)

        with h5py.File(self.filename, "a") as data:
            if self.method == "POD":
                self.runPOD(data, props, globdat)
            elif self.method == "LinearManifold":
                self.runLinMan(data, props, globdat)
            elif self.method == "QuadraticManifold":
                self.runQuadMan(data, props, globdat)
            else:
                raise NotImplementedError(f"{self.method} is not implemented")

        globdat.active = False

        self.writeFooter(globdat)

    def runPOD(self, data: h5py.File, props: Any, globdat: Any) -> None:
        """Construct a POD basis from the stored state snapshots.

        Args:
            data: Open HDF5 file containing at least the ``state`` dataset.
            props: Solver configuration properties. Not used directly.
            globdat: Global data object. Not used directly.

        Notes:
            The snapshot matrix is read from ``state`` and decomposed as
            ``X.T = U S V^T``. The left singular vectors are stored as ROM
            modes and the singular values are stored as ``eigenvals`` for later
            inspection.
        """
        snapshots = np.asarray(data["state"][:])
        sample_count, dof_count = snapshots.shape

        logger.info("    snapshot count   : %i", sample_count)
        logger.info("    state dimension  : %i", dof_count)

        modes, singular_values, _ = svd(snapshots.T)
        mode_count = modes.shape[1]

        logger.info("    mode count       : %i", mode_count)

        if singular_values.size > 0:
            logger.info("    leading sigma    : %6.4e", singular_values[0])
            if singular_values.size > 1:
                logger.info("    second sigma     : %6.4e", singular_values[1])

            total_energy = float(np.dot(singular_values, singular_values))
            if total_energy > 0.0:
                captured = np.cumsum(singular_values * singular_values) / total_energy
                logger.info("    energy(1 mode)   : %6.4f", captured[0])
                if singular_values.size > 1:
                    logger.info("    energy(2 modes)  : %6.4f", captured[1])

        if "modes" in data:
            del data["modes"]
        if "eigenvals" in data:
            del data["eigenvals"]

        data.create_dataset("modes", modes.shape, dtype="f", data=modes)
        data.create_dataset(
            "eigenvals",
            singular_values.shape,
            dtype="f",
            data=singular_values,
        )

        logger.info("    stored dataset   : modes %s", modes.shape)
        logger.info("    stored dataset   : eigenvals %s", singular_values.shape)

    def runLinMan(self, data: h5py.File, props: Any, globdat: Any) -> None:
        """Placeholder for linear-manifold ROM basis construction.

        Args:
            data: Open HDF5 file used for ROM data storage.
            props: Solver configuration properties.
            globdat: Global data object.

        Raises:
            NotImplementedError: Always, because the method is not implemented.
        """
        raise NotImplementedError("Linear Manifold Method is not yet implemented")

    def runQuadMan(self, data: h5py.File, props: Any, globdat: Any) -> None:
        """Placeholder for quadratic-manifold ROM basis construction.

        Args:
            data: Open HDF5 file used for ROM data storage.
            props: Solver configuration properties.
            globdat: Global data object.

        Raises:
            NotImplementedError: Always, because the method is not implemented.
        """
        raise NotImplementedError(
            "Quadratic Manifold Method is not yet implemented"
        )