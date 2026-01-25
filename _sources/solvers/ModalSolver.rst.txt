===========
ModalSolver
===========

The ``ModalSolver`` computes eigenmodes via the generalized eigenvalue problem
for the assembled stiffness and mass matrices. It stores eigenvectors and
values in the global data and terminates.

Overview
--------

Solver type: ``ModalSolver``

- Problem: solves ``K v = \lambda M v``
- Outputs: eigenvectors and eigenvalues
- Termination: deactivates model after solve

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
     - Must be set to ``"ModalSolver"``

Examples
--------

- Modal analysis in small plate or frame models

See Also
--------
- :doc:`../solvers`
- :doc:`DynEigSolver <DynEigSolver>`
