=====
Truss
=====

The ``Truss`` element models axial members with only axial stiffness and
no bending resistance, suitable for pin-jointed truss structures.

--------
Overview
--------

Element type: ``Truss``

Supports:
- Axial-only stiffness
- Pin-jointed structures and Riks analysis examples

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
     - Must be set to ``"Truss"``
   * - ``material``
     - Material block with axial properties (e.g., ``E``, area via section or mesh).

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``section``
     - Optional section/area definition.

--------
Examples
--------

Representative example using ``Truss``:
- ``examples/ch04/ShallowtrussRiks.pro``
