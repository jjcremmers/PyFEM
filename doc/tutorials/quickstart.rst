Quickstart
==========

In order to test whether everything is installed properly, the following two simulations can be run.

Simple example
--------------

In the directory ''examples/ch02'' the script ''PatchTest.py'' can be executed from a terminal (or DOS-shell) by typing:

  python PatchTest.py

In Windows, this script can also be executed by double-clicking the icon.

PyFEM example
-------------

The full finite element code PyFEM can be run by typing ''pyfem'' in the terminal. In directory ''examples/ch04'' for example,
the input file ''ShallowTrussRiks.pro'' is processed by typing::

  pyfem ShallowTrussRiks.pro

Here, ''ShallowTrussRiks.pro'' is the input file, which by definition ends with ''.pro''. When it is
opened in a text editor, it looks as follows::

  input = "ShallowTrussRiks.dat";

  TrussElem  = 
  {
    ....
  };

  SpringElem =
  {
    ....
  };

  solver =
  {
    ....
  };

  outputModules = ["graph"];

  graph = 
  {
    ....
  };
  

The dots indicate lines that have been omitted in this example.
The first argument in the ''.pro''-file specifies the input file, which contains the positions of 
the nodes, the element connectivity and the boundary conditions. The structure of this file, which normally has 
the extension ''.dat'', is as follows::

  <Nodes>
    1   0.0 0.0 ;
    2 -10.0 0.0 ;
    3  10.0 0.0 ;
    4   0.0 0.5 ;
  </Nodes>

  <Elements>
    1 'TrussElem' 2 4 ;
    2 'TrussElem' 3 4 ;
    3 'SpringElem' 1 4 ;
  </Elements>

  <NodeConstraints>
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;

    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 0.0;
  </NodeConstraints>

  <ExternalForces>
    v[4] = -100.0 ;
  </ExternalForces>

The nodes are defined between the labels ''<Nodes>'' and ''</Nodes>''. The first number 
indicates the node identification number. The remaining numbers denote the coordinates in x-, y-, and in the 
case of a three dimensional simulation, the z-direction. For example, node 2 has the x and y coordinates (-10,0).

The element connectivity is given after the tag ''<Elements>''. 
The first number indicates the element ID number. The string refers to the name of the
element model this element belongs to. The remaining numbers are the nodes that are used to construct 
the element. In this example, the first element is of the type ''TrussElem'' and is supported
by nodes 2 and 4.

The boundary conditions and applied loads are specified next. The node constraints are given after the label \texttt{<NodeConstraints>}. 
In this example, the displacement components ''u'' and
''v'' of nodes 1,2 and 3 have a prescribed value of 0.0. The external forces are specified in a 
similar manner in the field ''<ExternalForces>''. Here, a unit external force with magnitude -100.0 is
added to node number 4 in the direction that corresponds to the ''v'' displacement. Hence, this force is
acting in the negative y-direction.

The parameters of the finite element model are specified in the ''.pro'' file, 
in the fields ''TrussElem'' and ''SpringElem'', which refer to the labels used in the element connectivity
description::

  TrussElem  = 
  {
    type = "Truss";
    E    = 5e6;
    Area = 1.0;
  };

  SpringElem =
  {
    type = "Spring";
    k    = 100.0;
  };

The elements denoted by the label ''TrussElem'' are of the type ''Truss''. This model requires two
additional parameters, the Young's modulus of the material ''E'' and the area of the cross-section \texttt{Area}.
The label ''SpringElem'' denote elements of the type ''Spring''. Here, one additional parameter is required: the
spring stiffness ''k''. A detailed overview of the element types and the corresponding parameters can be found in Section
of this manual.

The parameters of the solver are defined next::

  solver =
  {
    type = 'RiksSolver';

    fixedStep = true;
    maxLam    = 10.0; 
  };

The solver is of the type ''RikSolver''. The two additional parameters specify that the magnitude of the path-parameter is
constant (''fixedStep = true'') and that the simulation is stopped when the load parameter $\lambda$ reaches a value of 10.0. 
A detailed overview of available solver types and their parameters is given in Section xx.

Finally, the results of the simulation can be stored and visualised in several ways. To this end, a chain of output modules can be 
specified. In this example, the results are stored in a load-displacement curve in the module ''GraphWriter''::

  outputModules = ["graph"];

  graph =
  {
    type     = "GraphWriter";
    onScreen = true;

    columns = [ "disp" , "load" ];

    disp =
    {
      type   = "state";
      node   = 4;
      dof    = 'v';
      factor = -1.0;
    };
  
    load =
    { 
      type = "fint";
      node = 4;
      dof  = 'v';
    };
  };

In this example, two colums are stored: ''disp'', the displacement (''state'') of node 4 in the vertical 
direction and ''load'', the corresponding internal force. The parameter ''onScreen = true'' is used
to show the load-displacement curve on the screen during the simulation. By default, the results will be stored in a file called
''ShallowTrussRiks.out''. A description of all available output 
modules can be found in Section X.
