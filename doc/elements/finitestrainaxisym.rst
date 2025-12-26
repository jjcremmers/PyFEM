=====================
FiniteStrainAxiSym
=====================

The ``FiniteStrainAxiSym`` element performs axisymmetric finite strain analysis
for rotationally symmetric structures undergoing large deformations.

--------
Overview
--------

Element type: ``FiniteStrainAxiSym``

Supports:
- Axisymmetric kinematics with finite strains
- Nonlinear analysis for large displacements

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
     - Must be set to ``"FiniteStrainAxiSym"``
   * - ``material``
     - Material block appropriate for axisymmetric finite strain analysis.

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``method``
     - Lagrangian formulation selector if supported (e.g., TL/UL).

--------
Examples
--------

Representative example using ``FiniteStrainAxiSym``:
- ``examples/elements/axisymmetric/axisymmetric_finitestrain4.pro``
