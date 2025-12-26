=============
SandwichCore
=============

The ``SandwichCore`` model provides orthotropic-like stiffness tailored for
sandwich cores with dominant through-thickness properties and reduced in-plane
stiffness.

Overview
--------

Material type: ``SandwichCore``

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
     - Must be set to ``"SandwichCore"``
   * - ``E3``
     - Through-thickness Young's modulus
   * - ``G``
     - Base shear modulus

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``G13`` / ``G23``
     - Shear moduli (defaults to ``G`` if not set)
   * - ``factor``
     - Reduction factor for in-plane stiffness (default 0.001)

Examples
--------

- Sandwich panels modeled with SLS elements

See Also
--------
- :doc:`../elements/sls`
