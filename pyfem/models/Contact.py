# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import numpy as np
from numpy import append,repeat
from pyfem.models.BaseModel import BaseModel

#===============================================================================
#
#===============================================================================

class Contact(BaseModel):
    """
    Contact model for enforcing contact constraints in finite element simulations.

    Attributes:
        flag (bool): Indicates if contact is active.
        dispDofs (list[str]): Displacement DOF names.
        centre (np.ndarray): Centre of the contact surface.
        direction (np.ndarray): Contact movement direction.
        radius (float): Contact surface radius.
        penalty (float): Penalty parameter for contact enforcement.
    """

    def __init__(self, props: object, globdat: object) -> None:
        """
        Initialize a Contact model instance with given properties and global data.

        Args:
            props: Model-specific properties (should include centre, radius, penalty, direction, type).
            globdat: Global data/state object.
        """
        BaseModel.__init__(self, props, globdat)
        self.flag = False
        self.dispDofs = ["u", "v"]
        self.centre = [1., 1.]
        self.direction = [0.0, 0.0]
        self.radius = 10.
        self.penalty = 1.e6
        if self.type == "sphere":
            self.dispDofs = ["u", "v", "w"]
        self.flag = True
        self.centre = np.array(props.centre)
        self.radius = props.radius
        self.penalty = props.penalty
        self.direction = np.array(props.direction)

#-------------------------------------------------------------------------------
#  checkContact   (with flag)
#-------------------------------------------------------------------------------


    def getTangentStiffness(self, props: object, globdat: object, mbuilder: object) -> None:
        """
        Compute and assemble the contact tangent stiffness matrix contribution.

        Args:
            props: Global properties object.
            globdat: Global data/state object.
            mbuilder: MatrixBuilder instance for assembly.
        """
        print("Contact getTangentStiffness")
        centre = self.centre + globdat.lam * self.direction
        for nodeID in list(globdat.nodes.keys()):
            crd = globdat.nodes.getNodeCoords(nodeID)
            idofs = globdat.dofs.getForTypes([nodeID], self.dispDofs)
            crd += globdat.state[idofs]
            ds = crd - centre
            dsnorm = np.linalg.norm(ds)
            overlap = self.radius - dsnorm
            if overlap > 0:
                normal = ds / dsnorm
                mbuilder.B[idofs] += -self.penalty * overlap * normal
                mat = self.penalty * np.outer(normal, normal)
                mbuilder.append(mat, idofs)