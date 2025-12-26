================
KirchhoffBeam
================

The ``KirchhoffBeam`` element models slender beams using the Kirchhoff (Euler)
beam theory, suitable for frames and slender structures where shear deformation
can be neglected.

--------
Overview
--------

Element type: ``KirchhoffBeam``

Supports:
- Frame analysis with bending and axial deformation
- Suitable for slender beams (no shear deformation)
- 2D frame applications demonstrated in examples

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
     - Must be set to ``"KirchhoffBeam"``
   * - ``material``
     - Material block with elastic properties (e.g., ``E``, ``rho``) and
       geometric properties provided via mesh/section data.

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

Representative examples using ``KirchhoffBeam``:
- ``examples/ch09/KirchhoffEuler.pro``
- ``examples/ch09/KirchhoffEuler_01.pro``
- ``examples/ch09/KirchhoffEuler_1.pro``
- ``examples/ch09/FrameKirchhoff.pro``
