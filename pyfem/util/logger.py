# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

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


def separator(symbol: str = "-", level: str = "info" ) -> None:
    """
    Log a horizontal separator line consisting of two leading spaces followed by
    81 copies of `symbol`.

    The separator is emitted using the root logger's INFO level so it respects
    the currently configured logging output format and level.

    Args:
        symbol: Single-character string used to draw the separator. Defaults to '-'.
    """
    if level == "debug":
        logging.getLogger().debug(symbol * 81)
    else:       
        logging.getLogger().info(symbol * 81)

def logVariable(description: str, value: Any, width: int = 35, level: str = "info") -> None:
    """
    Log a variable in a formatted way with aligned colons.
    
    Prints in the format: "[description] ............. : [value]"
    where the colon is always at the same position determined by `width`.
    
    Args:
        description: Description text for the variable (left side).
        value: Value to display (right side of colon).
        width: Total width for description + dots before colon. Defaults to 35.
        level: Logging level ("info" or "debug"). Defaults to "info".
    
    Examples:
        >>> logVariable("Number of nodes", 100)
        Number of nodes ............... : 100
        >>> logVariable("Young's modulus", 210e9, width=40)
        Young's modulus ........................ : 210000000000.0
    """
    # Calculate number of dots needed to reach the desired width
    num_dots = width - len(description)
    
    # Ensure at least one space and one dot
    if num_dots < 2:
        num_dots = 2
    
    # Create the formatted string
    formatted = f"{description} {'.' * num_dots} : {value}"
    
    # Log at the appropriate level
    if level == "debug":
        logging.getLogger().debug(formatted)
    else:
        logging.getLogger().info(formatted)


def logHeader(left_text: str, right_text: str, width: int = 35, level: str = "info") -> None:
    """
    Log two strings with aligned spacing similar to logVariable.
    
    Prints in the format: "[left_text] ............. : [right_text]"
    where the colon is always at the same position determined by `width`.
    
    Args:
        left_text: Text to display on the left side.
        right_text: Text to display on the right side of colon.
        width: Total width for left text + dots before colon. Defaults to 35.
        level: Logging level ("info" or "debug"). Defaults to "info".
    
    Examples:
        >>> logHeader("Parameter", "Value")
        Parameter ......................... : Value
        >>> logHeader("Solver type", "Nonlinear", width=40)
        Solver type ............................ : Nonlinear
    """
    # Calculate number of dots needed to reach the desired width
    numSpaces = width - len(left_text)
    
    # Ensure at least one space and one dot
    if numSpaces < 2:
        numSpaces = 2
    
    # Create the formatted string
    formatted = f"{left_text} {' ' * numSpaces} : {right_text}"
    
    # Log at the appropriate level
    if level == "debug":
        logging.getLogger().debug(formatted)
        separator(level="debug")
    else:
        logging.getLogger().info(formatted)        
        separator()
