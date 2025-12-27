===========
PlaneStress
===========

The ``PlaneStress`` material model implements linear elastic behavior in 2D
plane stress. Use for thin structures where out-of-plane stress is negligible.
Outputs stress components: ``S11``, ``S22``, ``S12``.

Overview
--------

Material type: ``PlaneStress``

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
     - Must be set to ``"PlaneStress"``
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

- ``examples/ch02/PatchTest4.pro``
- ``examples/plate/platetest.pro``
- ``examples/plate/plate_test_07.pro``
- ``examples/plate/platedyn.pro``

See Also
--------
- :doc:`../elements/plate`
- :doc:`../elements/smallstraincontinuum`
