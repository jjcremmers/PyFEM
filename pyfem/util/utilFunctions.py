# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from math import cos, sin
from typing import Union

import numpy as np


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


def sin_over_x(x_val: Union[int, float]) -> float:
    """
    Return sin(x) / x using a series expansion near zero.

    Args:
        x_val: Input value.

    Returns:
        The value of sin(x) / x.
    """

    if abs(x_val) < 1.0e-12:
        return 1.0 - x_val * x_val / 6.0

    return sin(x_val) / x_val


def one_minus_cos_over_x2(x_val: Union[int, float]) -> float:
    """
    Return (1 - cos(x)) / x^2 using a series expansion near zero.

    Args:
        x_val: Input value.

    Returns:
        The value of (1 - cos(x)) / x^2.
    """

    if abs(x_val) < 1.0e-12:
        return 0.5 - x_val * x_val / 24.0

    return (1.0 - cos(x_val)) / (x_val * x_val)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def array2tensor(a):

    b = np.zeros((3, 3))
    
    b[0, 0] =  0.0
    b[1, 0] =  a[2]
    b[2, 0] = -a[1]

    b[0, 1] = -a[2]
    b[1, 1] =  0.0
    b[2, 1] =  a[0]

    b[0, 2] =  a[1]
    b[1, 2] = -a[0]
    b[2, 2] =  0.0
    
    return b
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------  

def tensor2array(b):

    a = np.zeros(3)
    
    a[0] = b[2, 1]
    a[1] = b[0, 2]
    a[2] = b[1, 0]
    
    return a   
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def tosint(t):

    assert t >= 0.0, "Input t should be non-negative"
    
    if t <= 0.5:
        t2 = t * t
        result = 1.0 + t2 / 6.0
        t4 = t2 * t2
        result += 7.0 * t4 / 360.0
        t6 = t4 * t2
        result += 31.0 * t6 / 15120.0
        t8 = t4 * t4
        result += 127.0 * t8 / 604800.0
        t10 = t6 * t4
        result += 73.0 * t10 / 3421440.0
        t12 = t8 * t4
        result += 1414477.0 * t12 / 653837184000.0
        t14 = t8 * t6
        result += 8191.0 * t14 / 37362124800.0
        t16 = t8 * t8
        result += 16931177.0 * t16 / 762187345920000.0
        t18 = t10 * t8
        result += 5749691557.0 * t18 / 2554547108585472000.0
    else:
        result = t / np.sin(t)
    return result
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
def tantot(t):

    assert t >= 0.0, "Input t should be non-negative"
    
    if t <= 0.5:
        t2 = t * t
        result = 1.0 + t2 / 3.0
        t4 = t2 * t2
        result += 2.0 * t4 / 15.0
        t6 = t4 * t2
        result += 17.0 * t6 / 315.0
        t8 = t4 * t4
        result += 62.0 * t8 / 2835.0
        t10 = t6 * t4
        result += 1382.0 * t10 / 155925.0
        t12 = t8 * t4
        result += 21844.0 * t12 / 6081075.0
        t14 = t8 * t6
        result += 929569.0 * t14 / 638512875.0
        t16 = t8 * t8
        result += 6404582.0 * t16 / 10854718875.0
        t18 = t10 * t8
        result += 443861162.0 * t18 / 1856156927625.0
    else:
        result = np.tan(t) / t
    return result

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def sintot(t):

    assert t >= 0.0, "Input t should be non-negative"
    
    if t <= 0.5:
        t2 = t * t
        result = 1.0 - t2 / 6.0
        t4 = t2 * t2
        result += t4 / 120.0
        t6 = t4 * t2
        result -= t6 / 5040.0
        t8 = t4 * t4
        result += t8 / 362880.0
        t10 = t6 * t4
        result -= t10 / 39916800.0
        t12 = t8 * t4
        result += t12 / 6227020800.0
        t14 = t8 * t6
        result -= t14 / 1307674368000.0
        t16 = t8 * t8
        result += t16 / 355687428096000.0
    else:
        result = np.sin(t) / t
    return result

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def costot2(t):

    assert t >= 0.0, "Input t should be non-negative"
    
    if t <= 0.5:
        t2 = t * t
        result = 0.5 - t2 / 24.0
        t4 = t2 * t2
        result += t4 / 720.0
        t6 = t4 * t2
        result -= t6 / 40320.0
        t8 = t4 * t4
        result += t8 / 3628800.0
        t10 = t6 * t4
        result -= t10 / 479001600.0
        t12 = t8 * t4
        result += t12 / 87178291200.0
        t14 = t8 * t6
        result -= t14 / 20922789888000.0
        t16 = t8 * t8
        result += t16 / 6402373705728000.0
    else:
        result = (1 - np.cos(t)) / (t * t)
    return result

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def sintot3(t):

    assert t >= 0.0, "Input t should be non-negative"
    
    if t < 0.5:
        t2 = t * t
        result = 1.0 / 6.0 - t2 / 120.0
        t4 = t2 * t2
        result += t4 / 5040.0
        t6 = t4 * t2
        result -= t6 / 362880.0
        t8 = t4 * t4
        result += t8 / 39916800.0
        t10 = t6 * t4
        result -= t10 / 6227020800.0
        t12 = t8 * t4
        result += t12 / 1307674368000.0
        t14 = t8 * t6
        result -= t14 / 355687428096000.0
        t16 = t8 * t8
        result += t16 / 121645100408832000.0
    else:
        result = 1.0 / (t * t) - np.sin(t) / (t * t * t)
    return result
