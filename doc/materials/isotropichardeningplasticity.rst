=============================
IsotropicHardeningPlasticity
=============================

The ``IsotropicHardeningPlasticity`` model describes J2 plasticity with
isotropic hardening. Requires elastic properties and a yield stress, with a
hardening law defined via the `Hardening` utility.

Overview
--------

Material type: ``IsotropicHardeningPlasticity``

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
     - Must be set to ``"IsotropicHardeningPlasticity"``
   * - ``E``
     - Young's modulus
   * - ``nu``
     - Poisson's ratio
   * - ``syield``
     - Initial yield stress

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - Hardening law parameters
     - e.g., ``H``, tabulated functions, or other settings consumed by `Hardening`

Examples
--------

- ``examples/materials/plasticity/dogbone.pro``

See Also
--------
- :doc:`../elements/smallstraincontinuum`
- :doc:`../elements/finitestraincontinuum`
