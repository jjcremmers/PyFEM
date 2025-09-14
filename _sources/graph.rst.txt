========================
Graph output writer
========================

The output is stored in a multi column file by this writer. By defaults,
the file is stored with the extnsion `.out`. The first two columns of
the output can be shown on the screen as a curve during the simulation.
Please note that this option is only working when PyFEM is used in
a Linux environment.

-----------
Summary
-----------

+---------------------+-------------------------------------+
| Name                | GraphWriter                         |
+---------------------+-------------------------------------+
| Source              | pyfem/io/GraphWriter.py             |
+---------------------+-------------------------------------+

-----------------------
Mandatory parameters:
-----------------------
+---------------------+-------------------------------------+
| **Name**            | **Description**                     |
+---------------------+-------------------------------------+
| columns             | Array of strings indicating the     |
|                     | column that will be stored. For     |
|                     | each column, the type of data,      |
|                     | and if needed, the node, degree of  |
|                     | freedom and scaling factor needs    |
|                     | to be specified.                    |
+---------------------+-------------------------------------+
| type                | Type of data. This can be either    |
|                     | state, velo, fint, stress, etc.     |
+---------------------+-------------------------------------+
| node                |  Node ID.                           |
+---------------------+-------------------------------------+
| dof                 | Degree of freedom. This is most     |
|                     | likely ’u’ or ’v’                   |
+---------------------+-------------------------------------+

--------------------------
Optional parameters
--------------------------

+---------------------+-------------------------------------+
| **Name**            | **Description**                     |
+---------------------+-------------------------------------+
| factor              | The scaling factor for the output.  |
|                     | The default value is 1.0.           |
+---------------------+-------------------------------------+
| onScreen            | When set to true the first two      |
|                     | columns will be shown on the        |
|                     | screen. The default value is false. |
+---------------------+-------------------------------------+

--------------------------
Examples
--------------------------

The graph writer is used to plot the load-displacement curve in the
example in `examples/ch03/NewtonRaphson.pro'. In this case, the external load, 
the load factor `lambda` is plotted against the vertical displacement of node 21.

First the graph writer is called in the list of output modules ::

  outputModules = ["mesh" , "graph" ]; 

where `mesh` refers to the VTK output writer, which is discussed here, and 'graph' refers to this
`GraphWriter`. The following code block in the input-file describes the input paramers. ::

  graph =
  {
    type = "GraphWriter";
    onScreen = true;

    columns = [ "disp" , "load" ];

    disp = 
    {
      type = "state";
      node = 21;
      dof  = "v";
    };

    load = 
    {
      type = "lam";
    };
  };

In this example, the output file consists of two columns, labeled 'disp' and 'load'. The
first colum contains the displacement of node 21 in the y-direction. In the code, the displacement
is stored as 'state', whereas the 'v' degree-of-freedom is the displacement in the vertical direction.
In a similar fashion, the temperature or a rotation of a specific node can be plotted by
using the degree-of-freedom names, 'temp', 'rx', 'ry' or 'rz'.

The argument 'onScreen = true;' indicates that the results are plotted on screen during the simulation. Please not that this
only works in a Linux environment.




















