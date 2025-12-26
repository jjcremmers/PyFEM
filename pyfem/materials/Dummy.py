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

from typing import Tuple
from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import eye, dot, ndarray


class Dummy(BaseMaterial):
    """
    Dummy material model for interface elements.
    
    This is a simple linear elastic material model used primarily for testing
    and interface elements. It provides a diagonal stiffness matrix scaled by
    a material parameter D.
    
    Attributes
    ----------
    H : ndarray
        Material stiffness matrix (diagonal).
    outLabels : list of str
        Labels for output variables (normal and tangential tractions).
    D : float
        Material stiffness parameter (inherited from props).
    
    Notes
    -----
    The material supports 2D (rank=2) and 3D (rank=3) interface elements with
    normal and tangential traction components.
    """

    def __init__(self, props) -> None:
        """
        Initialize the Dummy material model.
        
        Parameters
        ----------
        props : object
            Properties object containing material parameters. Must have:
            - rank : int
                Spatial dimension (2 or 3).
            - D : float
                Stiffness parameter.
        """
        BaseMaterial.__init__(self, props)

        if props.rank == 2:
            self.H = self.D * eye(2)
            self.outLabels = ["Tn", "Ts"]
        elif props.rank == 3:
            self.H = self.D * eye(3)
            self.outLabels = ["Tn", "Ts1", "Ts2"]

    def getStress(self, deformation) -> Tuple[ndarray, ndarray]:
        """
        Compute stress (traction) and material tangent matrix.
        
        Parameters
        ----------
        deformation : object
            Deformation object containing the strain field.
            Must have a strain attribute (ndarray).
        
        Returns
        -------
        sigma : ndarray
            Stress (traction) vector.
        H : ndarray
            Material tangent stiffness matrix.
        """
        sigma = dot(self.H, deformation.strain)

        self.outData = sigma

        return sigma, self.H

    def getTangent(self) -> ndarray:
        """
        Get the material tangent stiffness matrix.
        
        Returns
        -------
        ndarray
            Material tangent stiffness matrix H.
        """
        return self.H

