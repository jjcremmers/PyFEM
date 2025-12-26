===========
RiksSolver
===========

The ``RiksSolver`` (arc-length method) follows equilibrium paths through
limit points by augmenting the Newton–Raphson iterations with an additional
constraint on arc-length. It adapts step size to balance convergence effort.

Overview
--------

Solver type: ``RiksSolver``

- Method: arc-length with predictor–corrector iterations
- Control: maintains a constraint on combined displacement and load increments
- Step size: adapts ``factor`` per cycle toward an optimal iteration count
- Termination: stops when ``lam > maxLam`` or a cycle cap is reached

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
     - Must be set to ``"RiksSolver"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``tol``
     - Convergence tolerance (default: ``1.0e-5``)
   * - ``optiter``
     - Target iterations per step to tune step factor (default: ``5``)
   * - ``iterMax``
     - Maximum iterations per cycle (default: ``10``)
   * - ``fixedStep``
     - Use fixed step size when ``true`` (default: ``false``)
   * - ``maxFactor``
     - Maximum cumulative factor to prevent step shrinkage (default: ``1.0e20``)
   * - ``maxLam``
     - Maximum load factor (default: ``1.0e20``)

Examples
--------

- ``examples/ch04/ShallowtrussRiks.pro``

See Also
--------
- :doc:`../solvers`
- Elements: :doc:`../elements/truss`, :doc:`../elements/smallstraincontinuum`
