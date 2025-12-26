=================
TimoshenkoBeam
=================

The ``TimoshenkoBeam`` element models beams including shear deformation effects,
providing improved accuracy for moderately thick beams compared to Kirchhoff
beam theory.

--------
Overview
--------

Element type: ``TimoshenkoBeam``

Supports:
- Bending with shear deformation
- Frame analysis where shear effects are non-negligible

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
     - Must be set to ``"TimoshenkoBeam"``
   * - ``material``
     - Material block with elastic properties (e.g., ``E``, ``G``, ``rho``).

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``section``
     - Optional section properties reference (if supported by input).

--------
Examples
--------

Representative examples using ``TimoshenkoBeam``:
- ``examples/ch09/FrameTimoshenko.pro``
