"""High-level API for programmatic control of a PyFEM analysis.

This module provides the `PyFEMAPI` class which encapsulates the complete
calculation machinery (input reader, solver, output manager) and exposes a
step-wise interface so callers can advance the analysis one step at a time or
drive the entire run to completion.

Usage:
    api = PyFEMAPI('mycase.pro')
    while api.is_active:
        api.step()
    results = api.get_results()
"""

from __future__ import annotations
from pathlib import Path
from typing import Any, Optional, Tuple, Union

from pyfem.io.InputReader import InputRead
from pyfem.io.InputReader import InputReader
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver import Solver


class PyFEMAPI:
    """Programmatic API to run PyFEM analyses step-by-step.

    The class can be constructed with either:
      - a path to a `.pro` file (string or Path), OR
      - a pre-built `(props, globdat)` tuple returned by `InputRead`.

    Methods:
        step(): advance one solver/output cycle
        run_all(): run until `globdat.active` becomes False
        is_active: boolean property exposing `globdat.active`
        get_results(): return a light-weight results dict (placeholder)
        close(): finalize and free resources
    """

    def __init__(self, props: Union[str, Path, Tuple[Any, Any]]) -> None:
        """Initialize API and internal components.

        Args:
            props: Either the path to a `.pro` file (or file basename), or a
                tuple `(props, globdat)` as returned by `InputRead`.
        """
        if isinstance(props, tuple) and len(props) == 2:
            self.props, self.globdat = props
        else:
            # Accept either a filename string or Path; InputRead expects the
            # pro filename as first argument.
            fname = str(props)
            self.props, self.globdat = InputRead(fname)

        # Instantiate solver and output manager
        self.solver = Solver(self.props, self.globdat)
        self.output = OutputManager(self.props, self.globdat)

    @property
    def isActive(self) -> bool:
        return bool(getattr(self.globdat, 'active', False))

    def step(self , nCyc: int = 1 ) -> None:
        """Perform a single solver step (by default) followed by output processing ."""

        for iCyc in range(nCyc):
            if not self.is_active:
                return
            self.solver.run(self.props, self.globdat)
            self.output.run(self.props, self.globdat)

    def runAll(self) -> None:
        """Run steps until the analysis completes."""
        while self.is_active:
            self.step()

    def getResults(self) -> Any:
        """Return a lightweight results container.

        This is a small convenience wrapper; consumers can directly inspect
        `self.globdat` for detailed state if needed.
        """
        return {
            'active': self.is_active,
            'globdat': self.globdat,
            'props': self.props,
        }

    def close(self) -> None:
        """Finalize the analysis and close global data resources."""
        try:
            self.globdat.close()
        except Exception:
            # Best-effort close; do not propagate to callers
            pass


def run(props: Union[str, Path, Tuple[Any, Any]]) -> Any:
    """Convenience function mirroring the old `run` behavior used in scripts.

    Creates a `PyFEMAPI` instance, runs the analysis to completion and returns
    the results dictionary.
    """
    api = PyFEMAPI(props)
    api.run_all()
    return api.get_results()

