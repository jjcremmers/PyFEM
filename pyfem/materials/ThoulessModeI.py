# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Any, Mapping, Tuple

import numpy as np

from pyfem.materials.BaseMaterial import BaseMaterial


class ThoulessModeI(BaseMaterial):
    """Mode-I cohesive zone model following Thouless.

    Parameters
    ----------
    props : Mapping[str, Any]
        Material properties expected by ``BaseMaterial`` to populate the
        attributes used here, notably ``Gc``, ``d1d3``, ``d2d3`` and ``Tult``.
    """

    def __init__(self, props: Mapping[str, Any]) -> None:
        """Initialize the material model and derived parameters.

        The base class processes ``props`` and defines attributes like ``Gc``,
        ``d1d3``, ``d2d3`` and ``Tult`` which are used to compute the internal
        parameters ``d1``, ``d2`` and ``d3``.
        """
        # Call the BaseMaterial constructor
        BaseMaterial.__init__(self, props)

        self.d3: float = 2.0 * self.Gc / ((-self.d1d3 + self.d2d3 + 1.0) * self.Tult)
        self.d1: float = self.d1d3 * self.d3
        self.d2: float = self.d2d3 * self.d3
        self.dummy: float = self.Tult / self.d1

        # Set the labels for the output data in this material model
        self.outLabels: list[str] = ["Tn", "Ts"]

    def getStress(self, deformation: Any) -> Tuple[np.ndarray, np.ndarray]:
        """Compute traction and consistent algorithmic tangent.

        Expects ``deformation.strain`` to provide the normal opening at index 0.
        Only the normal traction component is modeled in this material; the
        shear traction remains zero.

        Parameters
        ----------
        deformation : Any
            Object with attribute ``strain`` (NumPy array-like) containing the
            normal opening at index 0.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            A tuple ``(stress, tang)`` where ``stress`` is a length-2 array
            ``[Tn, Ts]`` and ``tang`` is the 2x2 algorithmic tangent. Only the
            ``(0, 0)`` component is non-zero in this model.

        Notes
        -----
        Also updates ``self.outData`` with the computed ``stress``.
        """
        stress = np.zeros(2)
        tang = np.zeros((2, 2))

        eps_n = deformation.strain[0]
        if eps_n < self.d1:
            stress[0] = self.dummy * eps_n
            tang[0, 0] = self.dummy
        elif self.d1 <= eps_n < self.d2:
            stress[0] = self.Tult
            tang[0, 0] = 0.0
        elif self.d2 <= eps_n < self.d3:
            stress[0] = self.Tult * (1.0 - (eps_n - self.d2) / (self.d3 - self.d2))
            tang[0, 0] = self.Tult * (-1.0) / (self.d3 - self.d2)
        else:
            stress[0] = 0.0
            tang[0, 0] = 0.0

        self.outData = stress

        return stress, tang
