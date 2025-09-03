=====
Plate
=====

The plate element is an implementation of the Kirchhoff-Love finite 
element as presented in Chapter 3 of the book XX. The element can 
be used to model flat structures in the ``x-y`` plane as 3, 4, 6 and 8 node 
elements with 5 degrees of freedom per node. Three translation degrees of 
freedom (``u`` , ``v`` and ``w``) which describe the displacements in the 
``x`` , ``y`` and ``z`` direction, respectively, and 2 rotational degrees 
of freedom (``rx`` and ``ry``) that represent the rotations around the 
``x`` and ``y`` axis, respectively.

---------
Arguments
---------

Element type: Plate



---------------------------
Example: Isotropic material
---------------------------

The plate element can be used to model a thin-walled structure, made of 
a single, isotropic material. In the following example, a plate with a 
thickness of 1.2 mm is considered, which is made of aluminium. The Young's
modulus is equal to :math:`E=72` GPa, the Poisson ratio :math:`\nu=0.3` and
:math:`\rho=2780` kg/m3.

The block in the input file that describes this plate is given below. Please 
not that all dimensions are in mm, kg and Pa.::

  PlateElem = 
  {
    type = "Plate";

    material = 
    {
      E    = 72e9;
      nu   = 0.3;
      rho  = 2780.;
    };

    thickness = 0.0012;
  };

This is the end.

---------------------------
Example: Layered composite
---------------------------

Alternatively, the element can be used to model flat, composite structures.
In the following example, a composite consisting of 5 layers is modeled, with
the following stacking sequence:

.. math::
  \lbrack 0_w , 0 , 90 , 0 , 0_w\rbrack
  
where :math:`0_w` is a woven layer thickness 0.22 mm with the following 
properties: :math:`E_1=10` GPa, :math:`E_2=10` GPa, :math:`\nu_{12}=0.25` and 
:math:`G_{12}=45` GPa. The three centre layers are made of a UD composite
with thickness 0.22 mm and properties: :math:`E_1=10` GPa, 
:math:`E_2=10` GPa, :math:`\nu_{12}=0.25` and :math:`G_{12}=45` GPa. 

These properties are given in the input file in the following way::
  
  PlateElem = 
  {
    type = "Plate";

    materials = [ "Woven" , "UD" ];

    layers = ["W" , "C0" , "C90" , "C0" , "W" ];

    Woven =
    {
      E1   = 1.e6;
      E2   = 0.5e5;
      nu12 = 0.3;
      G12  = 1.0e6;
      rho  = 1.0e3;
    };

    UD =
    {
      E1   = 1.e6;
      E2   = 0.5e5;
      nu12 = 0.3;
      G12  = 1.0e6;
      rho  = 1.0e3;
    };
  
    W = 
    {
      material  = "Woven";
      theta     = 0.;
      thickness = 0.05;
    };

    C0 =
    {
      material  = "UD";
      theta     = 0.;
      thickness = 0.05;   
    };
  
    C90 =
    {
      material  = "UD";
      theta     = 90.;
      thickness = 0.05;
    };
  };

This is the end.
