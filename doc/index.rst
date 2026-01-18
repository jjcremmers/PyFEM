.. PyFEM documentation master file, created by
   sphinx-quickstart on Mon Dec 30 18:39:21 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyFEM: A Python Finite Element Code
====================================

This is the user manual for PyFEM, a Python-based finite element code that
accompanies the book:

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
`Non-Linear Finite Element Analysis of Solids and
Structures <https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449>`__
John Wiley and Sons, 2012, ISBN 978-0470666449

.. figure:: https://media.wiley.com/product_data/coverImage300/47/04706664/0470666447.jpg
   :width: 200px
   :alt: The cover of Non-linear Finite element Analysis of Solids and Structures

The code is open source and intended for educational and scientific
purposes only. If you use PyFEM in your research, the developers would
be grateful if you could cite the book in your work.

Goals and Scope
---------------

PyFEM aims to provide a clear, well-documented reference implementation for
nonlinear finite element analysis, suited for teaching, prototyping, and
reproducible research. The code emphasizes readability over micro-optimizations
and includes a growing set of elements, material models, solvers, and I/O
modules to cover common solid mechanics problems.

How to Cite
-----------

If PyFEM contributes to a publication, please cite the textbook above and
reference PyFEM (with the commit/tag or release) to ensure reproducibility.
When applicable, include the specific modules (elements, materials, solvers)
used in your study.

J.J.C. Remmers (2026). PyFEM `<https://github.com/jjcremmers/PyFEM>`__

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel (2012)
`Non-Linear Finite Element Analysis of Solids and Structures <https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449>`__
John Wiley and Sons, 2012, ISBN 978-0470666449

License
-------

PyFEM is released under the MIT License, enabling broad use for education and research. 
For full license details, see the `LICENSE <../LICENSE>`__ file.

Under the MIT License terms, you may use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the software, provided that you
include the copyright notice and permission notice. The software is provided
"as is" without warranty of any kind.

Important
~~~~~~~~~

- The official PyFEM repository and its releases are licensed under the MIT
  License as described above (see `LICENSE <../LICENSE>`__).
- If you obtain PyFEM as part of a third-party distribution, fork, or bundled
  project, the licensing of that distribution may differ; always consult that
  project's license in addition to the original MIT-licensed PyFEM sources.
- Regardless of the surrounding project license, retain the file headers with
  attribution and disclaimer text present at the top of each source file.
- If you plan to adopt MIT licensing for your own distribution or fork, please
  keep the file-level notices intact and include an MIT license file in your
  distribution.

Documentation Contents
======================

.. toctree::
   :maxdepth: 1

   installation/overview
   quickstart
   usermanual
   tutorials/overview
   develop/overview
   api   
   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
