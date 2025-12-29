About the code
==============

This is the user manual for PyFEM. This python-based finite element code
accompanies the book:

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
`Non-Linear Finite Element Analysis of Solids and
Structures <https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449>`__
John Wiley and Sons, 2012, ISBN 978-0470666449

.. figure:: https://media.wiley.com/product_data/coverImage300/47/04706664/0470666447.jpg
   :width: 200px
   :alt: The cover of Non-linear Finite element Anaylsis of Solids and Structures

   alt text

The code is open source and intended for educational and scientific
purposes only. If you use PyFEM in your research, the developers would
be grateful if you could cite the book in your work.

[![Cite](https://img.shields.io/badge/Cite-How%20to%20cite-blue.svg)](#how-to-cite)

Goals and scope
---------------

PyFEM aims to provide a clear, well-documented reference implementation for
nonlinear finite element analysis, suited for teaching, prototyping, and
reproducible research. The code emphasizes readability over micro-optimizations
and includes a growing set of elements, material models, solvers, and I/O
modules to cover common solid mechanics problems.

How to cite
-----------

If PyFEM contributes to a publication, please cite the textbook above and
reference PyFEM (with the commit/tag or release) to ensure reproducibility.
When applicable, include the specific modules (elements, materials, solvers)
used in your study.

J.J.C. Remmers (2025). PyFEM [https://github.com/jjcremmers/PyFEM](https://github.com/jjcremmers/PyFEM)

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel (2012)
[`Non-Linear Finite Element Analysis of Solids and Structures'](https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449)
John Wiley and Sons, 2012, ISBN 978-0470666449

License
-------

PyFEM uses a project-level license together with file-level notices:

- Project-level license: see the file [LICENSE](LICENSE) (currently GNU GPLv3).
- File-level notices: each source file includes attribution and disclaimer
   text from the textbook and project history. Please retain these notices in
   redistributed or modified files.

MIT licensing option (documentation)
------------------------------------

To support broad educational and research use, PyFEM may be distributed under
the MIT License while keeping the file-level attribution/disclaimer headers.
Under the MIT terms you may use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the software, provided that you include the
copyright notice and permission notice. The software is provided “as is”
without warranty of any kind.

Important
~~~~~~~~~

- The effective legal terms for any given distribution are determined by the
   project-level license in the repository you use (see [LICENSE](LICENSE)).
- Regardless of the project-level license, retain the file headers with
   attribution and disclaimer text present at the top of each source file.
- If you plan to adopt MIT licensing for your own distribution or fork, please
   keep the file-level notices intact and include an MIT license file in your
   distribution.
