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
