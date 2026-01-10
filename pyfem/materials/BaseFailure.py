# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import List, Tuple, Any


class BaseFailure:
    """
    Base class for failure criteria in material models.
    
    This class serves as the foundation for implementing various failure 
    criteria (e.g., Tresca, von Mises, Mohr-Coulomb) used in material models.
    Properties are dynamically assigned from the input parameter list.
    
    Attributes:
        Properties are dynamically set based on the input props parameter.
    """

    def __init__(self, props: List[Tuple[str, Any]]) -> None:
        """
        Initialize the BaseFailure instance.
        
        Parameters
        ----------
        props : List[Tuple[str, Any]]
            A list of tuples containing property names and their corresponding 
            values. Each tuple should be (name, value) where name is a string 
            and value can be of any type.
        
        Examples
        --------
        >>> props = [('strength', 250.0), ('factor', 1.5)]
        >>> failure = BaseFailure(props)
        >>> failure.strength
        250.0
        """
        for name, val in props:
            setattr(self, name, val)  
