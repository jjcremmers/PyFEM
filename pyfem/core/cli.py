# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""Command-line entry point for running a PyFEM analysis.

This module provides the `main` function used by the console script to
initialize input/output managers and the selected solver and then execute the
analysis loop until completion.

The function intentionally keeps orchestration logic minimal: input parsing
and object factories are delegated to `InputReader`, `Solver` and
`OutputManager`.
"""

import sys
from argparse import ArgumentParser, Namespace

from pyfem import __version__
from pyfem.io.InputReader   import InputRead
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver   import Solver


def parse_arguments(argv: list[str] | None = None) -> Namespace:
    """Parse command-line arguments for the PyFEM executable."""
    parser = ArgumentParser(
        prog="pyfem",
        description="Run a PyFEM analysis using the specified input file.",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        help="Path to the PyFEM .pro input file.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_option",
        help="Path to the PyFEM .pro input file.",
    )
    parser.add_argument(
        "-d",
        "--dump",
        dest="dump_file",
        help="Restart from a PyFEM dump file.",
    )
    parser.add_argument(
        "-p",
        "--param",
        dest="parameters",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="Override an input parameter. May be supplied multiple times.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"PyFEM {__version__}",
    )

    args = parser.parse_args(argv)

    if args.input_file and args.input_option:
        parser.error("provide the input file either positionally or with -i/--input, not both")

    args.input_file = args.input_option or args.input_file

    if not args.input_file and not args.dump_file:
        parser.error("an input file or dump file is required")

    return args


def main(argv: list[str] | None = None) -> None:
    """
    Run a PyFEM analysis from the command line.

    Args:
        argv: Optional list of command-line arguments. If `None`, the program
            will use `sys.argv`.

    Notes:
        This function is an orchestration entry point; it does not perform any
        heavy computation itself. The created `props` and `globdat` objects are
        used to instantiate the `Solver` and `OutputManager` which perform the
        analysis and output duties respectively.
    """

    args = parse_arguments(argv)
    props, globdat = InputRead(args.input_file, args.dump_file, args.parameters)

    solver = Solver(props, globdat)
    output = OutputManager(props, globdat)

    while globdat.active:
        solver.run(props, globdat)
        output.run(props, globdat)

    globdat.close()


if __name__ == "__main__":
    main(sys.argv[1:])
