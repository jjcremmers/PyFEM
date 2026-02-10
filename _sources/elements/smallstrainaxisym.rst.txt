====================
SmallStrainAxiSym
====================

The ``SmallStrainAxiSym`` element performs axisymmetric small strain analysis
for rotationally symmetric structures.

--------
Overview
--------

Element type: ``SmallStrainAxiSym``

Supports:
- Axisymmetric kinematics with small strains
- Linear and mildly nonlinear analysis

----------
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
     - Must be set to ``"SmallStrainAxiSym"``
   * - ``material``
     - Material block appropriate for axisymmetric small strain analysis.

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``rho``
     - Density for dynamic problems if applicable.

--------
Examples
--------

Representative example using ``SmallStrainAxiSym``:
- ``examples/elements/axisymmetric/axisymmetric_smallstrain4.pro``
