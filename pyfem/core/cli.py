"""Command-line entry point for running a PyFEM analysis.

This module provides the `main` function used by the console script to
initialize input/output managers and the selected solver and then execute the
analysis loop until completion.

The function intentionally keeps orchestration logic minimal: input parsing
and object factories are delegated to `InputReader`, `Solver` and
`OutputManager`.
"""

import sys
from typing import Any

from pyfem.io.InputReader   import InputReader
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver   import Solver


def main(argv: list[str] | None = None) -> None:
    """
    Run a PyFEM analysis from the command line.

    Args:
        argv: Optional list of command-line arguments. If `None`, the program
            will use `sys.argv`. The actual parsing and interpretation of the
            arguments is performed by `InputReader` so `main` simply forwards
            the arguments to that component.

    Notes:
        This function is an orchestration entry point; it does not perform any
        heavy computation itself. The created `props` and `globdat` objects are
        used to instantiate the `Solver` and `OutputManager` which perform the
        analysis and output duties respectively.
    """

    props, globdat = InputReader(sys.argv) 

    solver = Solver(props, globdat)
    output = OutputManager(props, globdat)

    while globdat.active:
        solver.run(props, globdat)
        output.run(props, globdat)

    globdat.close()


if __name__ == "__main__":
    main(sys.argv[1:])