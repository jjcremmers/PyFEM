# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

from typing import Any, Callable
import sys

import h5py
import numpy as np
from numpy import zeros
from scipy import sparse

from pyfem.fem.Assembly import assembleExternalForce, assembleTangentStiffness
from pyfem.util.BaseModule import BaseModule
from pyfem.util.logger import getLogger

logger = getLogger()


class ReducedOrderSolver(BaseModule):
    """Solve the nonlinear equilibrium equations in a reduced subspace.

    The solver projects the tangent matrix and residual onto a precomputed ROM
    basis and performs Newton iterations on the reduced coordinates. The full
    state vector is reconstructed after every reduced solve step.
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the reduced-order nonlinear solver.

        Args:
            props: Solver configuration properties.
            globdat: Global data object containing the DOF space and solver
                status.

        Raises:
            RuntimeError: If the loaded basis size does not match the number of
                global degrees of freedom.
        """
        self.tol = 1.0e-3
        self.iterMax = 10
        self.maxCycle = sys.maxsize
        self.maxLam = 1.0e20
        self.dtime = 1.0
        self.loadFunc = "t"
        self.loadCases = []
        self.modeCount = 5

        BaseModule.__init__(self, props)

        if self.maxLam > 1.0e19 and self.maxCycle == sys.maxsize:
            self.maxCycle = 5

        globdat.lam = 0.0
        globdat.solverStatus.dtime = self.dtime

        self.loadfunc: Callable[[float], float] = eval(
            "lambda t : " + str(self.loadFunc)
        )

        if hasattr(self, "loadTable"):
            self.maxCycle = len(self.loadTable)
            load_table = zeros(self.maxCycle + 1)
            load_table[1:] = self.loadTable
            self.loadTable = load_table

        with h5py.File(self.modes, "r") as h5file:
            self.V = h5file["modes"][:, :self.modeCount]

        if self.V.shape[0] != len(globdat.dofs):
            raise RuntimeError("Number of DOFs does not coincide")

        logger.info("Starting reduced-order solver .........")
        logger.info("    modes file       : %s", self.modes)
        logger.info("    basis size       : %i x %i", self.V.shape[0], self.V.shape[1])
        logger.info("    requested modes  : %i", self.modeCount)
        logger.info("    tolerance        : %6.4e", self.tol)
        logger.info("    max iterations   : %i", self.iterMax)
        logger.info("    time step        : %6.4e", self.dtime)
        logger.info("    max cycles       : %i", self.maxCycle)
        logger.info("    max lambda       : %6.4e", self.maxLam)
        logger.info("    load function    : %s", self.loadFunc)
        if hasattr(self, "loadTable"):
            logger.info("    load table size  : %i", len(self.loadTable) - 1)
        if self.loadCases:
            logger.info("    extra load cases : %s", ", ".join(self.loadCases))

    def run(self, props: Any, globdat: Any) -> None:
        """Advance one nonlinear ROM load step.

        Args:
            props: Solver configuration properties.
            globdat: Global data object containing the current state, solver
                status, and assembled model data.

        Raises:
            RuntimeError: If the Newton-Raphson iterations do not converge
                within ``iterMax`` iterations.
        """
        self.writeHeader()

        stat = globdat.solverStatus
        stat.increaseStep()
        stat.iiter = 0

        dof_count = len(globdat.dofs)
        state = globdat.state
        state_increment = globdat.Dstate

        state_increment[:] = zeros(dof_count)

        logger.info("Nonlinear ROM solver ............")
        logger.info("    =============================================")
        logger.info("    Load step %i", globdat.solverStatus.cycle)
        logger.info("    =============================================")
        logger.info("    active dofs       : %i", dof_count)
        logger.info("    reduced modes     : %i", self.V.shape[1])
        logger.info("    Newton-Raphson   : L2-norm residual")

        self.setLoadAndConstraints(globdat)
        tangent, internal_force = assembleTangentStiffness(props, globdat)

        self.setLoadAndConstraints(globdat)
        external_force = assembleExternalForce(props, globdat)

        error = 1.0
        initial_residual = self.V.T @ (external_force - internal_force)

        logger.info("    load factor       : %6.4e", globdat.lam)
        logger.info("    load increment    : %6.4e", globdat.dlam)
        logger.info(
            "    initial residual  : %6.4e",
            np.linalg.norm(initial_residual),
        )

        while error > self.tol:
            stat.iiter += 1

            reduced_tangent = self.V.T @ sparse.csr_matrix.dot(tangent, self.V)
            reduced_residual = self.V.T @ (external_force - internal_force)

            reduced_increment = np.linalg.solve(reduced_tangent, reduced_residual)
            full_increment = self.V @ reduced_increment

            state_increment[:] += full_increment[:]
            state[:] += full_increment[:]

            tangent, internal_force = assembleTangentStiffness(props, globdat)

            error = reduced_increment @ reduced_increment

            logger.info(
                "    Iteration %4i   : %6.4e (|dq|=%6.4e, |du|=%6.4e)",
                stat.iiter,
                error,
                np.linalg.norm(reduced_increment),
                np.linalg.norm(full_increment),
            )

            if stat.iiter == self.iterMax:
                raise RuntimeError("Newton-Raphson iterations did not converge!")

        globdat.elements.commitHistory()

        state_increment[:] = zeros(len(globdat.dofs))
        globdat.fint = internal_force

        logger.info("    converged in      : %i iterations", stat.iiter)
        logger.info("    final load factor : %6.4e", globdat.lam)
        logger.info("    final residual    : %6.4e", error)

        if stat.cycle == self.maxCycle or globdat.lam > self.maxLam:
            globdat.active = False

        if not globdat.active:
            logger.info("    solver status     : stopping after this step")

        self.writeFooter(globdat)

    def setLoadAndConstraints(self, globdat: Any) -> None:
        """Update load factors and prescribed constraints for the current step.

        Args:
            globdat: Global data object containing solver status and degree of
                freedom information.
        """
        if hasattr(self, "loadTable"):
            cycle = globdat.solverStatus.cycle

            globdat.lam = self.loadTable[cycle]
            globdat.dlam = self.loadTable[cycle] - self.loadTable[cycle - 1]
            globdat.dofs.setConstrainFactor(globdat.dlam)
            globdat.solverStatus.lam = globdat.lam
            return

        globdat.lam = self.loadfunc(globdat.solverStatus.time)
        lam0 = self.loadfunc(
            globdat.solverStatus.time - globdat.solverStatus.dtime
        )

        globdat.dlam = globdat.lam - lam0
        globdat.dofs.setConstrainFactor(globdat.dlam)
        globdat.solverStatus.lam = globdat.lam

        logger.debug("  ---- main load -------------------------")
        logger.debug("    loadFactor       : %4.2f", globdat.lam)
        logger.debug("    incr. loadFactor : %4.2f", globdat.dlam)

        for load_case in self.loadCases:
            load_props = getattr(self.myProps, load_case)
            loadfunc = eval("lambda t : " + str(load_props.loadFunc))
            lam = loadfunc(globdat.solverStatus.time)
            lam0 = loadfunc(globdat.solverStatus.time - globdat.solverStatus.dtime)
            dlam = lam - lam0

            globdat.dofs.setConstrainFactor(dlam, load_props.nodeTable)

            logger.debug("  ---- %s ---------------------", load_case)
            logger.debug("    loadFactor       : %4.2f", lam)
            logger.debug("    incr. loadFactor : %4.2f", dlam)