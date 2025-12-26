==============
BuckEigSolver
==============

The ``BuckEigSolver`` performs a buckling eigenvalue analysis. It solves for
the prebuckling static state, assembles the updated stiffness, and computes
buckling eigenvalues/eigenvectors.

Overview
--------

Solver type: ``BuckEigSolver``

- Prebuckling: solves ``K0 a = f_ext`` for the static state
- Buckling: computes eigenpairs from the generalized problem using the initial
  and updated stiffness
- Outputs: ``eigenvals`` (critical loads) and ``eigenvecs`` (buckling modes)

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
     - Must be set to ``"BuckEigSolver"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``tol``
     - Convergence tolerance (default: ``1.0e-3``)
   * - ``iterMax``
     - Maximum iterations (if applicable)

Examples
--------

- Buckling studies on frames or plates with appropriate boundary conditions

See Also
--------
- :doc:`../solvers`
- :doc:`../io/MeshWriter`
