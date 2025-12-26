===
SLS
===

The ``SLS`` (Solid-Like Shell) element is an advanced shell element formulation 
that combines the computational efficiency of shell elements with the generality 
of solid elements. It is particularly well-suited for modeling laminated composite 
structures, layered materials, and sandwich panels where through-thickness behavior 
is important.

The element uses a solid-like kinematic description with displacement degrees of 
freedom (``u``, ``v``, ``w``) and can model multiple layers with different material 
properties and fiber orientations. It employs static condensation to eliminate 
internal degrees of freedom, maintaining computational efficiency while capturing 
complex through-thickness effects.

--------
Overview
--------

Element type: ``SLS``

The element implements:

- **Multi-layer capability**: Supports layered structures with different materials, 
  thicknesses, and fiber orientations per layer
- **Displacement-only formulation**: Uses only displacement DOFs (``u``, ``v``, ``w``) 
  without rotational degrees of freedom
- **Static condensation**: Internal degrees of freedom are condensed out to improve 
  computational efficiency
- **Transverse shear deformation**: Accounts for through-thickness shear effects
- **Geometric nonlinearity**: Can handle moderate geometric nonlinearity
- **Composite materials**: Fully compatible with anisotropic and transversely 
  isotropic material models

The element is ideal for modeling:

- Laminated composite plates and shells
- Sandwich structures with face sheets and core
- Multi-material layered systems
- Structures requiring accurate through-thickness stress predictions

----------
Parameters
----------

Mandatory Parameters
~~~~~~~~~~~~~~~~~~~~

.. list-table:: 
   :widths: 25 75
   :header-rows: 1

   * - Parameter
     - Description
   * - ``type``
     - Must be set to ``"SLS"``
   * - ``material``
     - Material block defining constitutive behavior. For single-layer elements, specify one material. For multi-layer laminates, use ``MultiMaterial`` type.
       
       Common material types:
       
       * ``"Isotropic"``: For isotropic materials (parameters: ``E``, ``nu``, ``rho``)
       * ``"TransverseIsotropic"``: For unidirectional composites (parameters: ``E1``, ``E2``, ``nu12``, ``G12``, ``rho``)
       * ``"Orthotropic"``: For fully orthotropic materials (parameters: ``E1``, ``E2``, ``E3``, ``nu12``, ``nu13``, ``nu23``, ``G12``, ``G13``, ``G23``, ``rho``)
       * ``"MultiMaterial"``: For laminates with different materials per layer

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table:: 
   :widths: 25 75
   :header-rows: 1

   * - Parameter
     - Description
   * - ``theta``
     - Fiber orientation angle in degrees for single-layer elements with anisotropic materials. Specifies the angle between the fiber direction and the element local x-axis. Default is 0.0 if not specified.
   * - ``layers``
     - List of layer identifiers for multi-layer laminates. Each layer must be defined as a separate block containing:
       
       * ``thickness``: Layer thickness
       * ``theta``: Fiber orientation angle for this layer (in degrees)
       * ``material``: Material name (when using ``MultiMaterial`` type)
       
       Example layer definition structure::
       
         layers = ["layer1", "layer2", "layer3"];
         
         layer1 = 
         {
           thickness = 0.5;
           theta     = 0.0;
           material  = "mat1";
         };

--------
Examples
--------

Example 1: Single-Layer Isotropic Shell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simple single-layer shell element with isotropic material::

  SLSElem =
  {
    type = "SLS";

    material = 
    {
      type = "Isotropic";
      E    = 1.e6;
      nu   = 0.0;
      rho  = 1.11e3;
    };
  };

This configuration is used in: ``examples/elements/sls/sls_cantilever01.pro``

Example 2: Unidirectional Composite with Fiber Orientation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A single-layer composite with fibers oriented at an angle::

  SLSElem =
  {
    type = "SLS";

    material = 
    {
      type = "TransverseIsotropic";
      
      E1   = 1.e6;
      E2   = 5.e5;
      nu12 = 0.25;
      G12  = 4.e5;
      rho  = 1.1e3;
    };
    
    theta = 45.0;     // Fiber angle in degrees
  };

Example 3: Multi-Layer Laminate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A three-layer composite laminate with different fiber orientations::

  SLSElem0 =
  {
    type = "SLS";

    material = 
    {
      type = "TransverseIsotropic";
      
      E1   = 1.e6;
      E2   = 5.e5;
      nu12 = 0.25;
      G12  = 4.e5;
      rho  = 1.1e3;
    };
    
    theta = 0.0;     // 0-degree layer
  };

  SLSElem1 =
  {
    type = "SLS";

    material = 
    {
      type = "TransverseIsotropic";
      
      E1   = 1.e6;
      E2   = 5.e5;
      nu12 = 0.25;
      G12  = 4.e5;
      rho  = 1.2e3;
    };
    
    theta = 90.0;    // 90-degree layer
  };

  SLSElem2 =
  {
    type = "SLS";

    material = 
    {
      type = "TransverseIsotropic";
      
      E1   = 1.e6;
      E2   = 5.e5;
      nu12 = 0.25;
      G12  = 4.e5;
      rho  = 1.3e3;
    };
    
    theta = 0.0;     // 0-degree layer
  };

This configuration (0/90/0 laminate) is used in: 
``examples/elements/sls/sls_cantilever02.pro``

------------------
Additional Examples
------------------

The ``SLS`` element is demonstrated in several examples in the 
``examples/elements/sls/`` directory:

- ``sls_cantilever01.pro``: Basic single-layer cantilever
- ``sls_cantilever02.pro``: Multi-layer laminate configuration
- ``sls_cantilever03.pro``: Advanced layered structure
- ``sls_cantilever04.pro``: Complex laminate example
- ``sls_cantilever_dyn.pro``: Dynamic analysis with SLS elements

---------------
Special Features
---------------

**Static Condensation**
  The SLS element uses static condensation to eliminate internal degrees of 
  freedom, making it computationally efficient while maintaining accuracy for 
  through-thickness behavior.

**Through-Thickness Integration**
  Multiple integration points through the thickness of each layer provide 
  accurate stress distributions and allow for nonlinear material behavior 
  variation across the thickness.

**Composite Modeling**
  Particularly well-suited for modeling composite laminates where each layer 
  can have different material properties, fiber orientations, and thicknesses.

See Also
--------

- :doc:`materials` - Material models for composites
- :doc:`plate` - Alternative plate element formulation
- :doc:`smallstraincontinuum` - Solid continuum elements
- :doc:`tutorial1` - Introduction to PyFEM input files
