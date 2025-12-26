################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since        #
#  then augmented and maintained by Joris J.C. Remmers.                        #
#  All rights reserved.                                                        #
#                                                                              #
#  A github repository, with the most up to date version of the code,          #
#  can be found here:                                                          #
#     https://github.com/jjcremmers/PyFEM/                                     #
#     https://pyfem.readthedocs.io/                                            #	
#                                                                              #
#  The original code can be downloaded from the web-site:                      #
#     http://www.wiley.com/go/deborst                                          #
#                                                                              #
#  The code is open source and intended for educational and scientific         #
#  purposes only. If you use PyFEM in your research, the developers would      #
#  be grateful if you could cite the book.                                     #    
#                                                                              #
#  Disclaimer:                                                                 #
#  The authors reserve all rights but do not guarantee that the code is        #
#  free from errors. Furthermore, the authors shall not be liable in any       #
#  event caused by the use of the program.                                     #
################################################################################

import logging
from typing import Any


def setLogger(props: Any) -> logging.Logger:
    """
    Create and configure the root logger for the current analysis.

    The function expects a configuration object which may provide
    a `logger.level` attribute (e.g. a namespace or simple object). The
    returned logger is the root logger configured with a single stream
    handler according to the requested verbosity.

    Args:
        props: Configuration object (may have a `logger` attribute with `.level`).

    Returns:
        logging.Logger: Configured root logger.

    Raises:
        NotImplementedError: If an unsupported level string is provided.
    """
    
    
    level = "normal"
  

    # Default to INFO level
    level = getattr(getattr(props, "logger", None), "level", "info")

    if level not in ["normal", "info", "debug", "critical", "warning", "silent"]:
        raise NotImplementedError(
            'Logger level should be "normal", "info", "debug", "critical", "silent" or "warning"'
        )

    logger = logging.getLogger()          # root logger
    logger.propagate = False              # avoid duplicates via parent loggers

    # Only add a handler if none exists
    if not logger.handlers:
        handler = logging.StreamHandler()

        # Configure formatter and level
        if level == "debug":
            formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
            logger.setLevel(logging.DEBUG)
        elif level in ["critical", "silent"]:
            formatter = logging.Formatter('  %(message)s')
            logger.setLevel(logging.CRITICAL)
        elif level == "warning":
            formatter = logging.Formatter('  %(message)s')
            logger.setLevel(logging.WARNING)
        else:
            formatter = logging.Formatter('  %(message)s')
            logger.setLevel(logging.INFO)

        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger
  
  
def getLogger() -> logging.Logger:
    """
    Return the active root logger instance.

    Returns:
        logging.Logger: The active root logger.
    """

    return logging.getLogger()


def separator(symbol: str = "-") -> None:
    """
    Log a horizontal separator line consisting of two leading spaces followed by
    81 copies of `symbol`.

    The separator is emitted using the root logger's INFO level so it respects
    the currently configured logging output format and level.

    Args:
        symbol: Single-character string used to draw the separator. Defaults to '-'.
    """
    logging.getLogger().info('  ' + (symbol * 81))
