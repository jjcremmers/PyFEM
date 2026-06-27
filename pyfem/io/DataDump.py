# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Any
import pickle
from pyfem.util.BaseModule import BaseModule


class DataDump(BaseModule):
    """Module for saving simulation state to disk using pickle serialization.
    
    Periodically saves the properties and global data to a pickle file,
    allowing simulation restart or post-processing.
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the DataDump module.
        
        Args:
            props: Properties dictionary containing module configuration.
            globdat: Global data object containing simulation state.
        """
        self.prefix = globdat.prefix
        self.extension = ".dump"
        self.lastOnly = False

        BaseModule.__init__(self, props)
        
        if not hasattr(props, "interval"):
            self.interval = 1
            
        if self.lastOnly:
            self.interval = 1

    def run(self, props: Any, globdat: Any) -> None:
        """Save simulation state to a pickle file.
        
        Writes the current properties and global data to a pickle file
        at intervals specified by the interval attribute.
        
        Args:
            props: Properties dictionary containing simulation parameters.
            globdat: Global data object containing current simulation state.
        """
        cycle = globdat.solverStatus.cycle
        
        if cycle % self.interval == 0:
            self.writeHeader()
            
            data = {}
            data["props"] = props
            data["globdat"] = globdat
            
            if self.lastOnly:
                name = str(self.prefix + self.extension)
            else:
                name = str(self.prefix + "_" + str(cycle) + self.extension)
            
            pickle.dump(data, open(name, "wb"))
