================
StaggeredSolver
================

``StaggeredSolver`` splits the solution into two subproblems handled by
configured sub-solvers with disjoint dof sets. It alternates solves for each
subset within a cycle and can perform Newton iterations for sub-solvers of type
``Nonlinear``.

Overview
--------

Solver type: ``StaggeredSolver``

- Sub-solvers: two configured blocks (``solver1``, ``solver2``) with names,
  types, and ``dofTypes`` defining their scope
- Load control: shared ``loadFunc`` and optional ``loadCases`` applied via
  sub-solver constrainers
- Termination: when reaching ``maxCycle``

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
     - Must be set to ``"StaggeredSolver"``
   * - ``solver1`` / ``solver2``
     - Sub-blocks describing each solver (see below)

Sub-Solver Block Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``name``
     - Identifier used in logs
   * - ``type``
     - ``"Nonlinear"`` or other recognized type
   * - ``dofTypes``
     - List of dof names handled by this sub-solver

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``tol``
     - Convergence tolerance for nonlinear sub-solvers (default: ``1.0e-3``)
   * - ``iterMax``
     - Max iterations per sub-solver (default: ``10``)
   * - ``dtime``
     - Time step for load function evaluation (default: ``1.0``)
   * - ``maxCycle``
     - Maximum cycles (default: unlimited)
   * - ``maxLam``
     - Maximum load factor (default: ``1.0e20``)
   * - ``loadFunc``
     - String expression in time (e.g., ``"t"``)
   * - ``loadCases``
     - List of named cases with sub-blocks defining ``loadFunc`` and ``nodeTable``

Examples
--------

- Staggered coupling between displacement and phase-field dofs (conceptual)

See Also
--------
- :doc:`../solvers`
- :doc:`NonlinearSolver <NonlinearSolver>`
