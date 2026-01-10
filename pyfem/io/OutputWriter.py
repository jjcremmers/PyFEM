# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Any, Optional
from pyfem.util.BaseModule import BaseModule
from pyfem.util.dataStructures import Properties
from numpy import ndarray, zeros


class OutputWriter(BaseModule):
    """Output writer module for writing nodal data to files.

    This module handles writing nodal output data to files during finite element
    analysis. It can write output to both files and screen, and inherits from
    BaseModule to provide configuration capabilities.

    Attributes:
        prefix (str): Prefix for output filenames (typically globdat.prefix + "_glob").
        extension (str): File extension for output files (default: ".out").
        onScreen (bool): Flag to enable/disable screen output (default: False).
        filename (str): Full filename for output (prefix + extension).
    """

    def __init__(self, props: Properties, globdat: Any) -> None:
        """Initialize the OutputWriter module.

        Sets up the output file prefix, extension, and screen output flag.
        The filename is constructed from the globdat prefix with "_glob" suffix.

        Args:
            props: Properties object containing configuration parameters.
            globdat: Global data object containing solver state and configuration.
                    Must have a 'prefix' attribute.
        """
        # Set default output file prefix (globdat prefix + "_glob")
        self.prefix = globdat.prefix + "_glob"
        
        # Set default file extension
        self.extension = ".out"
        
        # Disable screen output by default
        self.onScreen = False

        # Initialize base module with properties
        BaseModule.__init__(self, props)

        # Construct filename if not already set by BaseModule from props
        if not hasattr(self, "filename"):
            self.filename = self.prefix + self.extension

    # -------------------------------------------------------------------------
    # run - Execute the output writing operation
    # -------------------------------------------------------------------------

    def run(self, props: Properties, globdat: Any) -> None:
        """Execute the output writing operation.

        This method is called by the solver to write nodal output data. It writes
        a header with the current cycle number, optionally prints nodes to screen,
        and always writes nodes to the output file.

        Args:
            props: Properties object containing configuration parameters.
            globdat: Global data object containing solver state, nodes, and output data.
                    Must have:
                    - solverStatus.cycle: Current analysis cycle number
                    - printNodes method: Method to write nodal data
        
        Returns:
            None
        """
        # Write header with current cycle number
        self.writeHeader(globdat.solverStatus.cycle)

        # Optionally print nodes to screen
        if self.onScreen:
            globdat.printNodes()

        # Always write nodes to output file
        globdat.printNodes(self.filename)
