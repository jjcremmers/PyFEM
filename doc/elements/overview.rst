Elements
========

PyFEM has quite a number of element formulations. In due time, they will
be described here.

Continuum Elements
------------------

Solid continuum elements for bulk materials under small or finite strains.
Use these when modeling 2D/3D solids where displacements are the primary
unknowns and the response is governed by a continuum constitutive law.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   finitestrainaxisym.rst
   finitestraincontinuum.rst
   smallstraincontinuum.rst
   smallstrainaxisym.rst   

Shell and Plate Elements
------------------------

Beam, plate/shell, interface, and axisymmetric elements, plus utilities.
Choose these for slender members, layered shells, cohesive interfaces, or
specialized kinematics outside standard continuum formulations.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   plate.rst
   sls.rst

Beam and Rod Elements
---------------------

Beam, plate/shell, interface, and axisymmetric elements, plus utilities.
Choose these for slender members, layered shells, cohesive interfaces, or
specialized kinematics outside standard continuum formulations.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   beamnl.rst
   kirchhoffbeam.rst
   spring.rst
   timoshenkobeam.rst
   truss.rst

Other Elements
--------------

Beam, plate/shell, interface, and axisymmetric elements, plus utilities.
Choose these for slender members, layered shells, cohesive interfaces, or
specialized kinematics outside standard continuum formulations.

.. toctree::
   :maxdepth: 1
   :titlesonly:
   
   interface.rst
