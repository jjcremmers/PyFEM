=========
Isotropic
=========

The ``Isotropic`` material model implements linear elastic, isotropic behavior
for 3D solids. It can be used with continuum elements in both small strain and
finite strain analyses where the constitutive response is governed by Young's
modulus and Poisson's ratio. Density can be provided for dynamic problems.

--------
Overview
--------

Material type: ``Isotropic``

Supports:
- 3D linear elastic isotropic behavior
- Compatible with small strain and finite strain kinematics (via element choice)
- Stress output: ``S11``, ``S22``, ``S33``, ``S23``, ``S13``, ``S12``

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
     - Must be set to ``"Isotropic"``
   * - ``E``
     - Young's modulus (units consistent with the model setup)
   * - ``nu``
     - Poisson's ratio (typically between 0 and 0.5)

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``rho``
     - Mass density used in dynamic analyses (if applicable)
   * - ``incremental``
     - Optional flag for incremental stress update; default behavior computes
       stress from total strain

--------
Examples
--------

Representative examples using ``Isotropic`` materials:
- ``examples/ch02/PatchTest8_3D.pro`` (3D small-strain continuum)
- ``examples/ch03/NewtonRaphson3D.pro`` (3D finite-strain continuum)
- ``examples/elements/sls/sls_cantilever01.pro`` (SLS shell with isotropic faces)

See Also
--------
- :doc:`../elements/smallstraincontinuum`
- :doc:`../elements/finitestraincontinuum`
- :doc:`../elements/sls`
