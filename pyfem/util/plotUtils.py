# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import numpy as np
from typing import Sequence, Tuple


def plotCurve(output: Sequence[Tuple[float, float]]) -> None:
    """Plot a 2D curve using pylab.

    The function accepts any sequence of (x, y) pairs and plots them
    with red circle markers connected by lines. Import of plotting
    functions is local to the function to avoid enforcing a global
    dependency on matplotlib when the module is imported.

    Args:
        output: Sequence of (x, y) coordinate pairs to plot.

    Returns:
        None
    """

    from pylab import plot, show, xlabel, ylabel

    plot([x[0] for x in output], [x[1] for x in output], "r-o")
    show()


def plotTime(t: float) -> str:
    """Format a time duration into a human-readable string.

    Small durations are shown in scientific or fixed-point seconds,
    larger durations use minutes and hours as appropriate.

    Args:
        t: Time duration in seconds.

    Returns:
        A formatted string representing the elapsed time.
    """

    if t < 0.1:
        return f"{t:.1e} sec."
    if t < 60.0:
        return f"{t:.3f} sec."
    if t < 3600.0:
        minutes = int(t // 60)
        seconds = t % 60
        return f"{minutes} min. {seconds:.2f} sec."

    hours = int(t // 3600)
    minutes = int((t % 3600) // 60)
    seconds = t % 60
    return f"{hours} hrs. {minutes} min. {seconds:.2f} sec."


