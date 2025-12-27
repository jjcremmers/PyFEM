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


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


class BaseModel:
    """Base class for additional PyFEM models (Contact and RVE).

    This abstract base class defines the interface for additional model types in PyFEM.
    Models implement specificboundary conditions, or multi-physics
    coupling strategies. Derived classes should override the __init__ and run
    methods to provide model-specific functionality.

    Examples of derived models include:
    - Constraint models for prescribed displacements
    - Load models for external forces
    - Contact models for mechanical interaction
    - RVE models for periodic boundary conditions
    """

    def __init__ ( self, props , globdat ):
        """Initialize the model with properties and global data.

        Parameters
        ----------
        props : Properties
            Model-specific properties and configuration parameters read from
            the input file. These typically include model type, parameters,
            and behavior settings.
        globdat : GlobalData
            Global data structure containing:
            - nodes: Node coordinates and groups
            - elements: Element connectivity and properties
            - dofs: Degree-of-freedom space and constraints
            - state: Current displacement/state vector
            - fint: Internal force vector
            - fhat: External force vector

        Notes
        -----
        This base implementation does nothing. Derived classes should call
        super().__init__() and then perform model-specific initialization.
        """

        for name, val in props:
            setattr(self, name, val)
      

