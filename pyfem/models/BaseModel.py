# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers


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
      

