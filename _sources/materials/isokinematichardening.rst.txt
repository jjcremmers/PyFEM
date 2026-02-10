=============================
IsotropicKinematicHardening
=============================

The ``IsotropicKinematicHardening`` model describes J2 plasticity with
kinematic hardening (backstress evolution). Requires elastic properties,
yield stress, and a kinematic hardening modulus.

Overview
--------

Material type: ``IsotropicKinematicHardening``

Parameters
----------

Mandatory Parameters
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Must be set to ``"IsotropicKinematicHardening"``
   * - ``E``
     - Young's modulus
   * - ``nu``
     - Poisson's ratio
   * - ``syield``
     - Yield stress
   * - ``hard``
     - Kinematic hardening modulus

Examples
--------

- Kinematic hardening demonstrations in continuum elements

See Also
--------
- :doc:`../elements/smallstraincontinuum`
- :doc:`../elements/finitestraincontinuum`
