# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

"""
LARC03 failure criterion for unidirectional composite plies.

This module provides a Python implementation of the quasi-2D LARC03
criterion translated to the PyFEM failure-model interface. The implementation
supports separate matrix tension/compression and fiber tension/compression
failure mechanisms, including a fiber-misalignment-based compression mode.

The criterion is intended for stresses expressed in the material axes:

- ``stress[0]``: ``sigma11``
- ``stress[1]``: ``sigma22``
- ``stress[2]`` or ``stress[5]``: ``tau12``

Required properties
-------------------
Elastic constants:
- ``E1``
- ``E2``
- ``G12``
- ``nu12``

Strengths:
- ``Xt`` or ``XT`` or ``F1t``
- ``Xc`` or ``XC`` or ``F1c``
- ``Yt`` or ``YT`` or ``F2t``
- ``Yc`` or ``YC`` or ``F2c``
- ``S`` or ``SL`` or ``Fs``

Either provide in-situ strengths directly:
- ``YT_is``
- ``SL_is``

or provide the fracture-energy-based data used to derive them:
- ``GIc_L``
- ``GIIc_L``
- ``ply_thickness``

Optional properties:
- ``etaL``
- ``g``
- ``alpha0`` in radians
- ``alpha0_deg`` in degrees
"""

from math import atan, cos, pi, sin, sqrt, tan
from typing import Any

from numpy import ndarray

from pyfem.materials.BaseFailure import BaseFailure


def _positive_part(value: float) -> float:
    """Return the Macaulay bracket of ``value``.

    Parameters
    ----------
    value : float
        Scalar quantity that is clipped to zero when negative.

    Returns
    -------
    float
        ``value`` when it is positive and ``0.0`` otherwise.
    """

    return max(value, 0.0)


class Larc03(BaseFailure):
    """LARC03 composite failure criterion for quasi-2D stress states.

    The implementation follows the ply-level LARC03 formulation for
    unidirectional composites in the material frame. It combines separate
    checks for matrix tension, matrix compression, fiber tension, and fiber
    compression. The fiber-compression branch uses the misaligned-frame
    transformation described in the original formulation.

    Notes
    -----
    The model is loaded through the material block using ``failureType``. All
    failure properties therefore live in the same property set as the parent
    material definition.
    """

    def __init__(self, props) -> None:
        """Initialize the LARC03 criterion from material properties.

        Parameters
        ----------
        props : Properties
            Material-property container passed in from the PyFEM material
            manager. Required properties are documented in the module-level
            docstring. Strength aliases from legacy input decks are accepted.

        Raises
        ------
        ValueError
            If required elastic, strength, or in-situ data is missing.
        """

        BaseFailure.__init__(self, props)

        self.E1 = self._require_property("E1")
        self.E2 = self._require_property("E2")
        self.G12 = self._require_property("G12")
        self.nu12 = self._require_property("nu12")

        self.Xt = self._require_any_property("Xt", "XT", "F1t")
        self.Xc = self._require_any_property("Xc", "XC", "F1c")
        self.Yt = self._require_any_property("Yt", "YT", "F2t")
        self.Yc = self._require_any_property("Yc", "YC", "F2c")
        self.S = self._require_any_property("S", "SL", "Fs")

        self.GIc_L = self._get_any_property("GIc_L")
        self.GIIc_L = self._get_any_property("GIIc_L")
        self.ply_thickness = self._get_any_property("ply_thickness")

        self.YT_is = self._get_any_property("YT_is")
        self.SL_is = self._get_any_property("SL_is")

        self.etaL = self._get_any_property("etaL")
        self.g = self._get_any_property("g", default=-1.0)

        if hasattr(self, "alpha0"):
            self.alpha0 = self.alpha0
        else:
            alpha0_deg = self._get_any_property("alpha0_deg", default=53.0)
            self.alpha0 = alpha0_deg * pi / 180.0

        self.nu21 = self.nu12 * self.E2 / self.E1

        self._compute_in_situ_strengths()
        self._compute_derived_parameters()

    def check(self, stress: ndarray, deformation) -> float:
        """Evaluate the LARC03 failure index.

        Parameters
        ----------
        stress : ndarray
            Stress vector in material coordinates. Supported layouts are
            ``[sigma11, sigma22, tau12]`` for 2D and
            ``[sigma11, sigma22, sigma33, tau23, tau13, tau12]`` for 3D.
        deformation : object
            Kinematic state passed through the failure interface. It is not
            used directly by this criterion.

        Returns
        -------
        float
            Governing LARC03 failure index. Values greater than or equal to
            ``1.0`` indicate failure initiation.
        """

        sigma11 = stress[0]
        sigma22 = stress[1]
        tau12 = stress[2] if len(stress) == 3 else stress[5]

        matrix_fi = 0.0
        fiber_fi = 0.0
        extra_fi = 0.0

        if sigma22 >= 0.0:
            matrix_fi = self._compute_matrix_tension_index(sigma22, tau12)
        else:
            matrix_fi = self._compute_matrix_compression_index(
                sigma22, tau12, self.alpha0
            )

        if sigma11 >= 0.0:
            fiber_fi = sigma11 / self.Xt
        else:
            phi = self._compute_misalignment_angle(sigma11, sigma22, tau12)
            sigma11_m, sigma22_m, tau12_m = self._transform_to_misaligned_frame(
                sigma11, sigma22, tau12, phi
            )

            if sigma22_m < 0.0:
                fiber_fi = (abs(tau12_m) + self.etaL * sigma22_m) / self.SL_is
            else:
                fiber_fi = self._compute_matrix_tension_index(sigma22_m, tau12_m)

            extra_fi = self._compute_matrix_compression_index(
                sigma22_m, tau12_m, 0.0
            )

        return max(matrix_fi, fiber_fi, extra_fi)

    def _require_property(self, name: str) -> Any:
        """Return a mandatory property value.

        Parameters
        ----------
        name : str
            Name of the property that must be present on the instance.

        Returns
        -------
        Any
            The requested property value.

        Raises
        ------
        ValueError
            If the property is not available.
        """

        if not hasattr(self, name):
            raise ValueError(f"Required property '{name}' not found for Larc03")
        return getattr(self, name)

    def _get_any_property(self, *names: str, default: Any = None) -> Any:
        """Return the first available property from a list of aliases.

        Parameters
        ----------
        *names : str
            Candidate property names in order of preference.
        default : Any, optional
            Value returned when none of the aliases exist.

        Returns
        -------
        Any
            The first matching property value, or ``default``.
        """

        for name in names:
            if hasattr(self, name):
                return getattr(self, name)
        return default

    def _require_any_property(self, *names: str) -> Any:
        """Return the first available alias and raise when all are missing.

        Parameters
        ----------
        *names : str
            Candidate property names for the same physical quantity.

        Returns
        -------
        Any
            The resolved property value.

        Raises
        ------
        ValueError
            If none of the provided aliases exists.
        """

        value = self._get_any_property(*names)
        if value is None:
            aliases = ", ".join(names)
            raise ValueError(
                f"Required property not found for Larc03. Expected one of: {aliases}"
            )
        return value

    def _compute_matrix_compression_index(
        self, sigma22: float, tau12: float, alpha: float
    ) -> float:
        """Compute the matrix-compression failure index.

        Parameters
        ----------
        sigma22 : float
            Transverse normal stress in material direction 2.
        tau12 : float
            In-plane shear stress in the material frame.
        alpha : float
            Fracture-plane angle in radians.

        Returns
        -------
        float
            Matrix-compression contribution to the failure index.
        """

        cos_alpha = cos(alpha)
        sin_alpha = sin(alpha)

        tau_eff_t = _positive_part(
            -sigma22 * cos_alpha * (sin_alpha - self.etaT * cos_alpha)
        )
        tau_eff_l = _positive_part(
            cos_alpha * (abs(tau12) + self.etaL * sigma22 * cos_alpha)
        )

        return (tau_eff_t / self.ST) ** 2 + (tau_eff_l / self.SL_is) ** 2

    def _compute_matrix_tension_index(self, sigma22: float, tau12: float) -> float:
        """Compute the matrix-tension failure index.

        Parameters
        ----------
        sigma22 : float
            Transverse tensile stress in material direction 2.
        tau12 : float
            In-plane shear stress in the material frame.

        Returns
        -------
        float
            Matrix-tension contribution to the failure index.
        """

        y_ratio = sigma22 / self.YT_is
        s_ratio = tau12 / self.SL_is

        return (1.0 - self.g) * y_ratio + self.g * y_ratio * y_ratio + s_ratio * s_ratio

    def _compute_misalignment_angle(
        self, sigma11: float, sigma22: float, tau12: float
    ) -> float:
        """Compute the fiber-misalignment angle for fiber compression.

        Parameters
        ----------
        sigma11 : float
            Longitudinal stress in the fiber direction.
        sigma22 : float
            Transverse stress in the material frame.
        tau12 : float
            In-plane shear stress in the material frame.

        Returns
        -------
        float
            Misalignment angle in radians.

        Raises
        ------
        ValueError
            If the denominator of the closed-form expression becomes too small.
        """

        denom = self.G12 + sigma11 - sigma22

        if abs(denom) < 1.0e-14:
            raise ValueError(
                "Denominator in Larc03 fiber misalignment angle is too small."
            )

        return (abs(tau12) + (self.G12 - self.Xc) * self.phi_c) / denom

    def _transform_to_misaligned_frame(
        self, sigma11: float, sigma22: float, tau12: float, phi: float
    ) -> tuple[float, float, float]:
        """Rotate stresses to the misaligned fiber frame.

        Parameters
        ----------
        sigma11 : float
            Longitudinal stress in the original material frame.
        sigma22 : float
            Transverse stress in the original material frame.
        tau12 : float
            In-plane shear stress in the original material frame.
        phi : float
            Misalignment angle in radians.

        Returns
        -------
        tuple[float, float, float]
            Rotated stresses ``(sigma11_m, sigma22_m, tau12_m)``.
        """

        cphi = cos(phi)
        sphi = sin(phi)

        sigma11_m = cphi * cphi * sigma11 + sphi * sphi * sigma22 + 2.0 * sphi * cphi * tau12
        sigma22_m = sphi * sphi * sigma11 + cphi * cphi * sigma22 - 2.0 * sphi * cphi * tau12
        tau12_m = -sphi * cphi * sigma11 + sphi * cphi * sigma22 + (cphi * cphi - sphi * sphi) * tau12

        return sigma11_m, sigma22_m, tau12_m

    def _compute_lambda22_0(self) -> float:
        """Return the transverse compliance combination used by LARC03."""

        return 2.0 * (1.0 / self.E2 - (self.nu21 * self.nu21) / self.E1)

    def _compute_lambda44_0(self) -> float:
        """Return the in-plane shear compliance used by LARC03."""

        return 1.0 / self.G12

    def _compute_in_situ_strengths(self) -> None:
        """Compute in-situ strengths when only fracture data is provided.

        The method leaves user-supplied ``YT_is`` and ``SL_is`` unchanged. When
        they are absent, it derives them from the fracture energies and ply
        thickness using the standard LARC03 compliance relations.

        Raises
        ------
        ValueError
            If neither direct in-situ strengths nor the required fracture-data
            set is available, or if the derived compliance terms are not
            positive.
        """

        if self.YT_is is not None and self.SL_is is not None:
            return

        if None in (self.ply_thickness, self.GIc_L, self.GIIc_L):
            raise ValueError(
                "Larc03 requires either YT_is and SL_is, or ply_thickness together "
                "with GIc_L and GIIc_L."
            )

        lambda22 = self._compute_lambda22_0()
        lambda44 = self._compute_lambda44_0()

        if lambda22 <= 0.0 or lambda44 <= 0.0:
            raise ValueError("Derived compliance terms for Larc03 must be positive.")

        self.YT_is = sqrt(8.0 * self.GIc_L / (pi * self.ply_thickness * lambda22))
        self.SL_is = sqrt(8.0 * self.GIIc_L / (pi * self.ply_thickness * lambda44))

    def _compute_derived_parameters(self) -> None:
        """Compute secondary LARC03 parameters from the input data.

        This includes the mixed-mode weighting parameter ``g``, the in-situ
        friction parameter ``etaL`` when it is not supplied, the transverse
        shear strength term ``ST``, and the critical kinking angle ``phi_c``.

        Raises
        ------
        ValueError
            If the kinking-angle expression becomes singular.
        """

        if self.g < 0.0:
            if self.GIc_L is not None and self.GIIc_L is not None:
                self.g = self.GIc_L / self.GIIc_L
            else:
                self.g = 1.12 * (self._compute_lambda22_0() / self._compute_lambda44_0())
                self.g *= (self.YT_is / self.SL_is) ** 2

        if self.etaL is None:
            cos_alpha = cos(self.alpha0)
            self.etaL = -(self.SL_is * cos(2.0 * self.alpha0)) / (
                self.Yc * cos_alpha * cos_alpha
            )

        self.ST = self.Yc * cos(self.alpha0) * (
            sin(self.alpha0) + cos(self.alpha0) / tan(2.0 * self.alpha0)
        )
        self.etaT = -1.0 / tan(2.0 * self.alpha0)

        ratio = self.SL_is / self.Xc
        denom = 2.0 * (ratio + self.etaL)
        inside = max(1.0 - 4.0 * (ratio + self.etaL) * ratio, 0.0)

        if abs(denom) < 1.0e-14:
            raise ValueError("Denominator in Larc03 phi_c computation is too small.")

        self.phi_c = atan((1.0 - sqrt(inside)) / denom)