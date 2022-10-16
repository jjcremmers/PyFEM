========================
Graph output writer
========================

The output is stored in a multi column file by this writer. The 
first two columns can be shown on the screen as a curve during the simulation.

+---------------------+-------------------------------------+
| Name                | GraphWriter                         |
+---------------------+-------------------------------------+
| Source              | pyfem/io/GraphWriter.py             |
+---------------------+-------------------------------------+

-----------------------
Mandatory parameters:
-----------------------

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
Optional parameters:
--------------------------

+---------------------+-------------------------------------+
| factor              | The scaling factor for the output.  |
|                     | The default value is 1.0.           |
+---------------------+-------------------------------------+
| onScreen            | When set to true the first two      |
|                     | columns will be shown on the        |
|                     | screen. The default value is false. |
+---------------------+-------------------------------------+