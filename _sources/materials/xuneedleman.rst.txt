============
XuNeedleman
============

``XuNeedleman`` is a cohesive interface material model defining traction–
separation behavior with parameters governing peak traction and fracture energy.
Outputs traction components (rank-dependent).

Overview
--------

Material type: ``XuNeedleman``

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
     - Must be set to ``"XuNeedleman"``
   * - ``Gc``
     - Fracture energy
   * - ``Tult``
     - Ultimate traction

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``q`` / ``r``
     - Shape parameters for the traction–separation law

Examples
--------

- ``examples/elements/interface/PeelTest.pro``
- ``examples/ch13/TractionOscillations.pro``

See Also
--------
- :doc:`../elements/interface`
