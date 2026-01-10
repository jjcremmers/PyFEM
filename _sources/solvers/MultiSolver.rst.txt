===========
MultiSolver
===========

``MultiSolver`` runs a sequence of solver blocks in order. When one solver
finishes (deactivates the model), ``MultiSolver`` activates the next solver and
continues until all are complete.

Overview
--------

Solver type: ``MultiSolver``

- Orchestration: imports and constructs solver classes by name
- Progression: advances to the next solver when ``globdat.active`` becomes
  false

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
     - Must be set to ``"MultiSolver"``
   * - ``solvers``
     - List of solver block names to run in order

Solver Blocks
-------------

For each name in ``solvers``, define a block with at least a ``type``:

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Solver class name under ``pyfem.solvers`` (e.g., ``NonlinearSolver``)
   * - Other parameters
     - Solver-specific options (see their documentation)

Examples
--------

- Run a nonlinear solve followed by a modal analysis:

.. code-block:: text

   solver =
   {
     type    = "MultiSolver";
     solvers = [ "nl", "modal" ];
   };

   nl = { type = "NonlinearSolver"; maxCycle = 10; };
   modal = { type = "DynEigSolver"; eigenCount = 5; };

See Also
--------
- :doc:`../solvers`
- Individual solver docs for block parameters
