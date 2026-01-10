# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import importlib
from typing import Any

from pyfem.util.dataStructures import Properties, GlobalData


class Solver:
    """
    Loader and base wrapper for solver implementations.

    This class dynamically imports the solver module specified in
    props.solver.type, instantiates the solver class and delegates the
    run call to that instance.

    Attributes:
        solver: Instantiated solver object (implementation-specific)
    """

    def __init__(self, props: Properties, globdat: GlobalData) -> None:
        """
        Load and instantiate the solver specified in properties.

        Args:
            props: Analysis properties (must contain a 'solver' attribute
                   with a 'type' field indicating the solver class name)
            globdat: GlobalData instance passed to the solver constructor

        Raises:
            ImportError: If the solver module or class cannot be found.
        """
        solverProps = getattr(props, "solver")
        solverType: str = solverProps.type

        try:
            mod = importlib.import_module(f"pyfem.solvers.{solverType}")
            SolverClass = getattr(mod, solverType)
        except ModuleNotFoundError as e:
            raise ImportError(
                f"Solver module 'pyfem.solvers.{solverType}' not found."
            ) from e
        except AttributeError as e:
            raise ImportError(
                f"Class '{solverType}' not found in module "
                f"'pyfem.solvers.{solverType}'."
            ) from e

        props.currentModule = "solver"
        self.solver: Any = SolverClass(props, globdat)

    def run(self, props: Properties, globdat: GlobalData) -> None:
        """
        Delegate execution to the selected solver instance.

        Args:
            props: Analysis properties
            globdat: GlobalData instance (may be modified in-place)
        """
        self.solver.run(props, globdat)
