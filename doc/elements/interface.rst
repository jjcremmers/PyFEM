================
Interface
================

The ``Interface`` element models cohesive or contact interfaces between continuum
parts. It is used to simulate delamination, peeling, and traction-separation
behaviors using appropriate interface constitutive laws.

--------
Overview
--------

Element type: ``Interface``

Supports:
- Traction–separation laws via material models (e.g., ``XuNeedleman``)
- 2D and 3D interface formulations depending on mesh and setup
- Nonlinear analysis for softening and debonding

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
     - Must be set to ``"Interface"``
   * - ``material``
     - Interface material model block defining cohesive law and parameters.

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``thickness``
     - Effective interface thickness (if required by material model).
   * - ``orientation``
     - Optional local orientation of the interface.

--------
Examples
--------

Representative examples using ``Interface``:
- ``examples/elements/interface/PeelTest.pro``
- ``examples/ch13/PeelTest.pro``
- ``examples/ch13/PeelTest60.pro``
- ``examples/ch13/PeelTest60pres.pro``
- ``examples/ch13/PeelTest80.pro``
- ``examples/ch13/PeelTest100.pro``
- ``examples/ch13/TractionOscillations.pro``
- ``examples/solver/dissipatedEnergySolver/delam_buckling100.pro``
- ``examples/solver/dissipatedEnergySolver/delam_buckling200.pro``

See Also
--------
- :doc:`materials` for cohesive laws (e.g., Xu–Needleman)
- :doc:`solvers` for nonlinear analysis options
