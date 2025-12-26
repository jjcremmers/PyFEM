======
BeamNL
======

The ``BeamNL`` element models beams with geometric nonlinearity, suitable for
buckling, large-displacement frame analysis, and dynamic eigenvalue problems.

--------
Overview
--------

Element type: ``BeamNL``

Supports:
- Nonlinear frame and column problems
- Buckling and dynamic analyses

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
     - Must be set to ``"BeamNL"``
   * - ``material``
     - Material block with elastic properties.

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``section``
     - Optional section properties reference.

--------
Examples
--------

Representative examples using ``BeamNL``:
- ``examples/ch09/Frame.pro``
- ``examples/elements/beam/verticalColumn.pro``
- ``examples/solver/buckEigSolver/eulerBuck.pro``
- ``examples/solver/dynEigSolver/beamDyn.pro``
