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

import numpy as np

def plotCurve( output: np.ndarray ) -> None:

    """
    Plots a curve based on the given output data points.

    Args:
        output (List[Tuple[float, float]]): A list of (x, y) data points to plot.

    Returns:
        None
    """

    from pylab import plot, show, xlabel, ylabel

    plot( [x[0] for x in output], [x[1] for x in output], 'r-o' )

    show()
  
          
def plotTime(t: float) -> str:
    
    """
    Formats a given time duration into a human-readable string.

    Args:
        t (float): Time duration in seconds.

    Returns:
        str: A formatted string representing the time in seconds, minutes, or hours.
    """
    if t < 0.1:
        return f"{t:.1e} sec."
    elif t < 60.0:
        return f"{t:.3f} sec."
    elif t < 3600.0:
        minutes = int(t // 60)
        seconds = t % 60
        return f"{minutes} min. {seconds:.2f} sec."
    else:
        hours = int(t // 3600)
        minutes = int((t % 3600) // 60)
        seconds = t % 60
        return f"{hours} hrs. {minutes} min. {seconds:.2f} sec."
        

