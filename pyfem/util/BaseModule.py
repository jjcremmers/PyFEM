# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""
Base module for PyFEM.

This module provides the BaseModule class, which serves as the base class for
modules in the PyFEM framework. It handles initialization of module properties
and provides utility methods for logging and timing module execution.

The BaseModule automatically detects whether it's a solver module and provides
appropriate logging levels. It also manages nested property hierarchies from
configuration objects.
"""

import time
from typing import Any, Optional
from pyfem.util.logger import getLogger,separator
from pyfem.util.plotUtils import plotTime

logger = getLogger()


class BaseModule:
    """
    Base class for all PyFEM modules.

    This class provides common functionality for modules including:
    - Automatic detection of module type (solver vs. non-solver)
    - Property management from configuration objects
    - Consistent logging with timing information
    - Support for hierarchical property organization

    Attributes:
        isSolver (bool): True if this module is a solver, False otherwise.
        type (str): The class name of the module.
        myProps: The properties dictionary for this module from the configuration.

    Args:
        props: A properties object containing configuration. Can have nested
            attributes or a currentModule attribute specifying the module name.
    """

    def __init__(self, props: Any) -> None:
        """
        Initialize the base module with properties from configuration.

        Attempts to load module-specific properties from the props object,
        supporting both single-level and nested property hierarchies. Also
        automatically detects if the module is a solver.

        Args:
            props (Any): Configuration object containing module properties.
                Can have a currentModule attribute specifying the module name,
                or will use the class name. Supports nested modules like
                "parent.child" in addition to flat module names.
        """
        self.isSolver = False

        if "Solver" in self.__class__.__name__.lower():
            self.isSolver = True

        if hasattr(props, 'currentModule') and hasattr(props, props.currentModule):
            currentModule = props.currentModule

            if 'solver' in currentModule.lower():
                self.isSolver = True
        elif hasattr(props, 'currentModule'):
            currentModule = props.currentModule

            if 'solver' in currentModule.lower():
                self.isSolver = True

        else:
            currentModule = self.__class__.__name__

            print(currentModule)
            if currentModule.endswith("olver") == "olver":
                currentModule = "solver"

                self.isSolver = True

        c = currentModule.split('.')

        if len(c) == 1:
            if hasattr(props, currentModule):
                self.myProps = getattr(props, currentModule)

                for name, val in self.myProps:
                    setattr(self, name, val)
        elif len(c) == 2:
            if hasattr(props, c[0]):
                p2 = getattr(props, c[0])
                if hasattr(p2, c[1]):
                    self.myProps = getattr(p2, c[1])

                    for name, val in self.myProps:
                        setattr(self, name, val)

        self.type = self.__class__.__name__

    def writeHeader(self, cycle: Optional[int] = None) -> None:
        """
        Write a header message at the start of module execution.

        Logs module name and step information. Solver modules use info-level
        logging with decorative separators, while non-solver modules use
        debug-level logging.

        Args:
            cycle (int, optional): The step/cycle number to include in the header.
                If provided, will be formatted as "step: {cycle}". Defaults to None.

        Side Effects:
            - Records start time in self.t0 for later timing calculations
            - Logs to logger.info or logger.debug depending on module type
        """
        self.t0 = time.time()
        cycleString = ""

        if cycle is not None:
            cycleString = "  step: " + str(cycle)

        if self.isSolver:
            separator()
            separator("=")
            logger.info("  " + self.type + cycleString)
            separator("=")
        else:
            separator(level = "debug" )
            logger.debug("  Module " + self.type)
            separator(level = "debug" )

    def writeFooter(self, globdat: Any) -> None:
        """
        Write a footer message at the end of module execution.

        Logs timing information for both the current step and total execution time.
        Solver modules use info-level logging, while non-solver modules use
        debug-level logging.

        Args:
            globdat (Any): Global data object containing startTime attribute
                for the overall simulation start time.

        Side Effects:
            - Logs step and total elapsed times via logger
            - Uses self.t0 (set by writeHeader) as start time reference
        """
        t1 = time.time()

        if self.isSolver:
            logger.info("    Elapsed time (this step).. : " + plotTime(t1 - self.t0))
            logger.info("    Total elapsed time........ : " + plotTime(t1 - globdat.startTime))
            separator()
        else:
            logger.debug("    Elapsed time (this step).. : " + plotTime(t1 - self.t0))
            logger.debug("    Total elapsed time........ : " + plotTime(t1 - globdat.startTime))      
            separator(level = "debug" )   

