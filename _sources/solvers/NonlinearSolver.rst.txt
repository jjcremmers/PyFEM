================
NonlinearSolver
================

The ``NonlinearSolver`` advances the analysis in cycles using a Newton–Raphson
procedure. At each cycle it assembles the tangent stiffness and internal force,
solves for the displacement increment, updates the state, and checks
convergence against a residual norm. Loads can be defined via a time-based
function or a tabulated sequence and may include multiple cases.

Overview
--------

Solver type: ``NonlinearSolver``

- Method: standard Newton–Raphson with residual norm check
- Termination: when convergence is reached each cycle; run stops at
  ``maxCycle`` or when ``lam > maxLam``
- Load control: ``loadFunc`` evaluated over time or ``loadTable`` per cycle;
  additional ``loadCases`` supported

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
     - Must be set to ``"NonlinearSolver"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``tol``
     - Convergence tolerance on residual norm (default: ``1.0e-3``)
   * - ``iterMax``
     - Maximum Newton iterations per cycle (default: ``10``)
   * - ``dtime``
     - Time increment per cycle used for evaluating ``loadFunc`` (default: ``1.0``)
   * - ``maxCycle``
     - Maximum number of cycles (default: derived; set explicitly to limit runs)
   * - ``maxLam``
     - Maximum load factor; run stops when exceeded (default: ``1.0e20``)
   * - ``loadFunc``
     - String expression in ``t`` (time), e.g., ``"t"`` or ``"sin(t)"``
   * - ``loadTable``
     - List of load factors per cycle; overrides ``loadFunc`` when present
   * - ``loadCases``
     - List of named sub-cases configured as properties on the solver block

Per-Load-Case Sub-Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using ``loadCases``, define a sub-block for each case name:

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``loadFunc``
     - String expression for the case’s load factor vs time
   * - ``nodeTable``
     - Node indices or group used to scope constraint factors for the case

Algorithm Notes
---------------

- Assembles ``K`` and ``f_int``; computes external force ``f_ext``.
- Solves ``K da = f_ext - f_int`` and updates ``a`` and ``Da``.
- Residual norm: if ``||f_ext|| < 1e-16``, uses ``||f_ext - f_int||``;
  otherwise uses normalized residual ``||f_ext - f_int|| / ||f_ext||``.
- On convergence: commits element history and writes output via configured
  I/O modules.

Examples
--------

Minimal ``.pro`` solver configuration:

.. code-block:: text

   solver =
   {
     type     = "NonlinearSolver";
     tol      = 1.0e-4;
     iterMax  = 20;
     dtime    = 1.0;
     maxCycle = 20;
     loadFunc = "t";
   };

Examples in the repository:

- ``examples/ch03/NewtonRaphson.pro``
- ``examples/ch03/NewtonRaphson3D.pro``

See Also
--------
- :doc:`../solvers`
- Elements: :doc:`../elements/smallstraincontinuum`, :doc:`../elements/finitestraincontinuum`
