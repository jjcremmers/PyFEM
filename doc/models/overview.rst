======
Models
======

PyFEM includes specialized computational models that extend beyond standard
finite element formulations. These models implement additional physics,
constraints, or multi-scale approaches for specific analysis types.

Overview
--------

Models in PyFEM provide:

- **Periodic boundary conditions** for representative volume elements (RVE)
- **Contact mechanics** for interface interactions
- **Multi-scale coupling** between different length scales
- **Special constraint systems** for specific applications

Models are configured in the ``.pro`` input file and run alongside the solver
to apply constraints, manage coupling, or enforce special conditions.

Configuration
-------------

Models are specified in the input file using a model block:

.. code-block:: text

   model_name = 
   {
     type = "ModelType";
     # model-specific parameters
   };

Multiple models can be active simultaneously, each handling different aspects
of the analysis.

Available Models
----------------

.. toctree::
   :maxdepth: 1
   :titlesonly:

   rve

See Also
--------

- :doc:`../elements/overview` - Element formulations
- :doc:`../materials/overview` - Material models
- :doc:`../solvers/overview` - Solution algorithms
