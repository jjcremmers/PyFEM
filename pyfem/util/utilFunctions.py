# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

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
