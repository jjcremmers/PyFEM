================
Developer's Guide
================

This guide provides comprehensive documentation for extending PyFEM with new
finite element formulations, material models, solution algorithms, and I/O
modules. Each section includes theoretical background, implementation examples,
and best practices.

Overview
--------

PyFEM is designed to be extensible, allowing researchers and developers to
implement new capabilities while leveraging existing infrastructure. The
modular architecture separates concerns into distinct components:

- **Elements**: Define finite element formulations and local computations
- **Materials**: Implement constitutive laws relating stress to strain
- **Solvers**: Orchestrate the solution process and convergence
- **I/O Modules**: Handle input/output and data management

All implementations follow object-oriented design patterns and integrate
seamlessly with PyFEM's assembly routines, DOF management, and output system.

Theoretical Foundation
----------------------

The implementations in PyFEM are based on the textbook:

.. note::
   **"Non-Linear Finite Element Analysis of Solids and Structures"**
   
   by R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
   
   John Wiley & Sons, 2012, ISBN 978-0470666449

Throughout the developer guides, specific equations and algorithms from this
book are referenced to provide theoretical context for implementations.

Development Topics
------------------

Element Development
~~~~~~~~~~~~~~~~~~~

Learn how to implement new finite element formulations including shape
functions, strain-displacement relations, and element matrix assembly.

:doc:`elements_dev`

Key topics:

- Element class structure and inheritance
- B-matrix computation for various element types
- Integration with material models
- Mass matrices for dynamic analysis
- Patch tests and validation

Material Model Development
~~~~~~~~~~~~~~~~~~~~~~~~~~

Implement constitutive laws defining stress-strain relationships for linear
elastic, elastoplastic, and cohesive zone models.

:doc:`materials_dev`

Key topics:

- Material class structure
- Elastic materials (plane stress, plane strain)
- Plasticity with hardening (von Mises, Drucker-Prager)
- Cohesive zone models for fracture
- State variable management
- Consistent tangent operators

Solver Development
~~~~~~~~~~~~~~~~~~

Create solution algorithms for static, dynamic, and stability analysis
including Newton-Raphson, arc-length, and explicit time integration.

:doc:`solvers_dev`

Key topics:

- Solver class structure and flow control
- Newton-Raphson iteration
- Arc-length methods for snap-back problems
- Explicit dynamics solvers
- Convergence criteria and error handling
- Load stepping strategies

I/O Module Development
~~~~~~~~~~~~~~~~~~~~~~

Develop input/output modules for reading problem definitions, writing results
in various formats (text, binary, VTK), and managing data flow.

:doc:`io_dev`

Key topics:

- I/O module structure and integration
- Text-based output writers
- Binary formats (HDF5)
- VTK files for visualization
- Time history tracking
- State saving and restoration

Getting Started
---------------

Prerequisites
~~~~~~~~~~~~~

Before developing new modules, familiarize yourself with:

1. **Python programming**: Object-oriented design, NumPy arrays
2. **Finite element theory**: Chapters 1-4 of the textbook
3. **PyFEM structure**: Run existing examples to understand workflow
4. **Version control**: Git basics for contributing code

Development Workflow
~~~~~~~~~~~~~~~~~~~~

1. **Study existing implementations**: Review similar elements/materials/solvers
2. **Design your module**: Plan class structure and methods
3. **Implement incrementally**: Start with simple cases, add complexity
4. **Test thoroughly**: Unit tests, patch tests, benchmark problems
5. **Document clearly**: Docstrings, equations, usage examples
6. **Contribute back**: Submit pull requests to the main repository

Code Style and Conventions
---------------------------

Naming Conventions
~~~~~~~~~~~~~~~~~~

- **Classes**: ``PascalCase`` (e.g., ``SmallStrainContinuum``)
- **Methods**: ``camelCase`` (e.g., ``getTangentStiffness``)
- **Variables**: ``camelCase`` or ``snake_case`` consistently
- **Constants**: ``UPPER_CASE``

Documentation
~~~~~~~~~~~~~

All classes and methods should include docstrings:

.. code-block:: python

   class MyElement(Element):
       """Brief description of element.
       
       Detailed explanation including theoretical basis,
       assumptions, and references to the textbook.
       
       Parameters:
           param1: Description
           param2: Description
       """
       
       def myMethod(self, arg1, arg2):
           """Brief method description.
           
           Args:
               arg1: Description
               arg2: Description
               
           Returns:
               Description of return value
           """
           pass

Type Hints
~~~~~~~~~~

Use type hints for clarity:

.. code-block:: python

   from typing import Tuple
   import numpy as np
   
   def getStress(self, strain: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
       """Compute stress and tangent."""
       return sigma, tangent

Testing Guidelines
------------------

All new modules should include:

1. **Unit tests**: Test individual methods in isolation
2. **Integration tests**: Test with complete analyses
3. **Validation tests**: Compare against analytical solutions
4. **Regression tests**: Ensure changes don't break existing functionality

Place tests in the ``test/`` directory:

.. code-block:: text

   test/
       test_my_element.py
       test_my_material.py
       test_my_solver.py

Contributing
------------

To contribute your developments to PyFEM:

1. Fork the repository on GitHub
2. Create a feature branch
3. Implement and test your module
4. Document thoroughly
5. Submit a pull request

See the PyFEM repository for detailed contribution guidelines:
https://github.com/jjcremmers/PyFEM

Additional Resources
--------------------

Example Implementations
~~~~~~~~~~~~~~~~~~~~~~~

Study these well-documented implementations as templates:

- **Element**: ``pyfem/elements/SmallStrainContinuum.py``
- **Material**: ``pyfem/materials/PlaneStress.py``
- **Solver**: ``pyfem/solvers/NonlinearSolver.py``
- **I/O Module**: ``pyfem/io/GraphWriter.py``

Utility Modules
~~~~~~~~~~~~~~~

Leverage existing utilities:

- ``pyfem/util/shapeFunctions.py``: Shape function evaluation
- ``pyfem/util/kinematics.py``: Kinematic transformations
- ``pyfem/fem/Assembly.py``: Matrix assembly routines
- ``pyfem/util/logger.py``: Logging functionality

Developer Documentation
~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   elements_dev
   materials_dev
   solvers_dev
   io_dev

References
----------

**Primary Reference**

*"Non-Linear Finite Element Analysis of Solids and Structures"*
by R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
John Wiley & Sons, 2012, ISBN 978-0470666449

**Additional References**

- Hughes, T.J.R. (2000). *The Finite Element Method*. Dover Publications.
- Belytschko, T., Liu, W.K., Moran, B. (2000). *Nonlinear Finite Elements for 
  Continua and Structures*. Wiley.
- Zienkiewicz, O.C., Taylor, R.L. (2000). *The Finite Element Method* (5th ed.). 
  Butterworth-Heinemann.

Support and Contact
-------------------

For questions about PyFEM development:

- GitHub Issues: https://github.com/jjcremmers/PyFEM/issues
- Email: Contact the PyFEM developers through the GitHub repository

When asking questions, please:

1. Provide a minimal example demonstrating the issue
2. Include PyFEM version and Python environment details
3. Reference specific sections of the developer guides
4. Show what you've already tried
