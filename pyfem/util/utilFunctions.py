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

from typing import Union


def macauley(x: Union[int, float]) -> Union[int, float]:
    """
    Compute the Macaulay bracket operation on a value.
    
    The Macaulay operation returns the input value when positive and zero
    when negative. Commonly used in mechanics for stress and strain calculations.
    
    Args:
        x: The input value (int or float)
        
    Returns:
        The Macauley value: x if x >= 0, else 0
        
    Examples:
        >>> macauley(5.0)
        5.0
        >>> macauley(-3.0)
        0
        >>> macauley(0.0)
        0.0
    """
    
    if x >= 0.:
        return x
    else:
        return 0


def sign(x: Union[int, float]) -> float:
    """
    Return the sign of a value.
    
    Returns 1.0 for positive or zero values, and -1.0 for negative values.
    
    Args:
        x: The input value (int or float)
        
    Returns:
        The sign: 1.0 if x >= 0, else -1.0
        
    Examples:
        >>> sign(5.0)
        1.0
        >>> sign(-3.0)
        -1.0
        >>> sign(0.0)
        1.0
    """
    
    if x < 0.:
        return -1.0
    else: 
        return 1.0
