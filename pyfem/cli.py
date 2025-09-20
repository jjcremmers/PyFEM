# pyfem/cli.py
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

from pyfem.io.InputReader   import InputReader
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver   import Solver


def main(argv: list[str] | None = None) -> None:

    props,globdat = InputReader( sys.argv )

    solver = Solver        ( props , globdat )
    output = OutputManager ( props , globdat )

    while globdat.active:
        solver.run( props , globdat )
        output.run( props , globdat )

    globdat.close()
    

if __name__ == "__main__":
    main(sys.argv[1:])

