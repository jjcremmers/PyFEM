======
Spring
======

The ``Spring`` element models axial spring behavior between nodes and can be
used to add elastic supports or connections in structural models.

--------
Overview
--------

Element type: ``Spring``

Supports:
- Axial spring stiffness
- Used alongside truss elements in Riks-type analyses

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
     - Must be set to ``"Spring"``
   * - ``material``
     - Spring stiffness and properties (depending on input format).

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``direction``
     - Optional directional definition if applicable.

--------
Examples
--------

Representative example using ``Spring``:
- ``examples/ch04/ShallowtrussRiks.pro``
