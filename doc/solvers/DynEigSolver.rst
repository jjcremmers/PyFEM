=============
DynEigSolver
=============

``DynEigSolver`` computes a selected number of eigenfrequencies and modes from
``K`` and ``M`` and reports the frequencies (rad/s and Hz). Eigenvectors are
stored in the global data for postprocessing.

Overview
--------

Solver type: ``DynEigSolver``

- Problem: solves ``K v = \lambda M v`` and reports ``\omega = sqrt(\lambda)``
- Outputs: ``eigenvecs`` and ``eigenvals`` (rad/s)
- Reporting: logs mode number, angular frequency, and frequency in Hz

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
     - Must be set to ``"DynEigSolver"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``eigenCount``
     - Number of modes to compute (default: ``5``)
   * - ``tol``
     - Convergence tolerance (if used by eigensolver)

Examples
--------

- Dynamic modal analysis using plate or continuum elements

See Also
--------
- :doc:`../solvers`
- :doc:`ModalSolver <ModalSolver>`
