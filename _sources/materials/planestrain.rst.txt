===========
PlaneStrain
===========

The ``PlaneStrain`` material model implements linear elastic behavior in 2D
plane strain. Use for thick structures or long domains where out-of-plane
strain is negligible. Outputs: ``S11``, ``S22``, ``S12``.

Overview
--------

Material type: ``PlaneStrain``

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
     - Must be set to ``"PlaneStrain"``
   * - ``E``
     - Young's modulus
   * - ``nu``
     - Poisson's ratio

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``rho``
     - Mass density for dynamic analyses

Examples
--------

- ``examples/ch05/StressWave_20x20.pro``
- ``examples/elements/interface/PeelTest.pro``
- ``examples/gmsh/two_fibres.pro``

See Also
--------
- :doc:`../elements/smallstraincontinuum`
- :doc:`../elements/finitestraincontinuum`
