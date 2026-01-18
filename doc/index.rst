.. PyFEM documentation master file, created by
   sphinx-quickstart on Mon Dec 30 18:39:21 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyFEM: A Python Finite Element Code
====================================

This is the user manual for PyFEM, a Python-based finite element code that
accompanies the book:

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel  
`Non-Linear Finite Element Analysis of Solids and Structures <https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449>`__  
John Wiley and Sons, 2012, ISBN 978-0470666449


.. raw:: html

  <div style="display: flex; align-items: center; justify-content: center; gap: 100px; margin-bottom: 1em;">
     <a href="https://media.wiley.com/product_data/coverImage300/47/04706664/0470666447.jpg" target="_blank">
       <img src="https://media.wiley.com/product_data/coverImage300/47/04706664/0470666447.jpg" alt="The cover of Non-linear Finite element Analysis of Solids and Structures" style="width:200px; box-shadow: 2px 2px 8px rgba(0,0,0,0.18); display:inline-block; vertical-align:middle;"/>
     </a>
     <a href="_static/book_page.png" target="_blank">
       <img src="_static/book_page.png" alt="Sample book page - Click to enlarge" style="width:200px; box-shadow: 2px 2px 8px rgba(0,0,0,0.18); display:inline-block; vertical-align:middle;"/>
     </a>
   </div>

.. raw:: html

   <p>
     <a href="https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449" target="_blank">
       <img alt="Wiley" src="https://img.shields.io/badge/Wiley-Buy-blue?logo=wiley" style="margin-right:5px;"/>
     </a>
     <a href="https://www.amazon.com/dp/0470666447" target="_blank">
       <img alt="Amazon" src="https://img.shields.io/badge/Amazon-Buy-orange?logo=amazon"/>
     </a>
     <a href="https://www.bol.com/nl/nl/p/nonlinear-finite-element-analysis-of-solids-and-structures/9200000002230150/?s2a=&bltgh=tyDCX4LwsfJHhGbeEVtF3A.4_25_26.27.FeatureOptionButton#productTitle" target="_blank">
       <img alt="Bol.com" src="https://img.shields.io/badge/Bol.com-Buy-blue?logo=bolcom"/>
     </a>
     <a href="https://www.barnesandnoble.com/w/non-linear-finite-element-analysis-of-solids-and-structures-ren-de-borst/1101190392" target="_blank">
       <img alt="Barnes & Noble" src="https://img.shields.io/badge/Barnes%20%26%20Noble-Buy-green?logo=barnesandnoble"/>
     </a>
     <a href="#" onclick="var bib=document.getElementById('bibtex-cite-book1'); bib.style.display=(bib.style.display==='block')?'none':'block'; return false;">
       <img alt="BibTeX" src="https://img.shields.io/badge/BibTeX-Cite-red?logo=latex"/>
     </a>
     <div id="bibtex-cite-book" style="display:none; background:#f8f8f8; border:1px solid #ddd; border-radius:6px; padding:12px; margin:10px 0; font-family:monospace; font-size:0.95em; position:relative;">
       <pre id="bibtex-text-book" style="margin:0; background:none; border:none; font-family:inherit; font-size:inherit; overflow-x:auto; white-space:pre;">
         @book{deBorst2012,
           title={Non-Linear Finite Element Analysis of Solids and Structures},
           author={de Borst, R. and Crisfield, M.A. and Remmers, J.J.C. and Verhoosel, C.V.},
           year={2012},
           publisher={John Wiley and Sons},
           isbn={978-0470666449}
         }
       </pre>
     </div>
   </p>

The code is open source and intended for educational and scientific
purposes only. If you use PyFEM in your research, the developers would
be grateful if you could cite the book in your work.

Goals and Scope
---------------

PyFEM aims to provide a clear, well-documented reference implementation for
nonlinear finite element analysis, suited for teaching, prototyping, and
reproducible research. The code emphasizes readability over micro-optimizations
and includes a growing set of :doc:`elements/overview`, 
:doc:`material models <materials/overview>`, :doc:`solvers/overview`, and 
:doc:`I/O modules <io/overview>` to cover common solid mechanics problems.

How to Cite
-----------

If PyFEM contributes to a publication, please cite the textbook above and
reference PyFEM (with the commit/tag or release) to ensure reproducibility.
When applicable, include the specific modules (elements, materials, solvers)
used in your study.


J.J.C. Remmers (2026). PyFEM `<https://github.com/jjcremmers/PyFEM>`__

.. raw:: html

   <p>
     <a href="#" onclick="var bib=document.getElementById('bibtex-cite-code'); bib.style.display=(bib.style.display==='block')?'none':'block'; return false;">
       <img alt="BibTeX" src="https://img.shields.io/badge/BibTeX-Cite-red?logo=latex"/>
     </a>
     <div id="bibtex-cite-book" style="display:none; background:#f8f8f8; border:1px solid #ddd; border-radius:6px; padding:12px; margin:10px 0; font-family:monospace; font-size:0.95em; position:relative;">
       <pre id="bibtex-text-book" style="margin:0; background:none; border:none; font-family:inherit; font-size:inherit; overflow-x:auto; white-space:pre;">
         @misc{PyFEM2026,
           author    = {J.J.C. Remmers},
           title     = {PyFEM},
           year      = {2026},
           howpublished = {\url{https://github.com/jjcremmers/PyFEM}},
           note      = {Python Finite Element Method code.}
         }
       </pre>
     </div>
   </p>

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel (2012)
`Non-Linear Finite Element Analysis of Solids and Structures <https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449>`__
John Wiley and Sons, 2012, ISBN 978-0470666449

.. raw:: html

   <p>
     <a href="#" onclick="var bib=document.getElementById('bibtex-cite-book2'); bib.style.display=(bib.style.display==='block')?'none':'block'; return false;">
       <img alt="BibTeX" src="https://img.shields.io/badge/BibTeX-Cite-red?logo=latex"/>
     </a>
     <div id="bibtex-cite-book" style="display:none; background:#f8f8f8; border:1px solid #ddd; border-radius:6px; padding:12px; margin:10px 0; font-family:monospace; font-size:0.95em; position:relative;">
       <pre id="bibtex-text-book" style="margin:0; background:none; border:none; font-family:inherit; font-size:inherit; overflow-x:auto; white-space:pre;">
         @book{deBorst2012,
           title={Non-Linear Finite Element Analysis of Solids and Structures},
           author={de Borst, R. and Crisfield, M.A. and Remmers, J.J.C. and Verhoosel, C.V.},
           year={2012},
           publisher={John Wiley and Sons},
           isbn={978-0470666449}
         }
       </pre>
     </div>
   </p>

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
