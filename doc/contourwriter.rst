++++++++++++++
Contour writer
++++++++++++++

The output of a selected number of nodes is stored in a multi column 
file by this writer. For each step in the simulation a new file will be 
written the following format:

``filename-contour-stepnumber.out``

The columns contain the nodeIDs, the x and y coordinates of the node, the
displacements, stresses and the tractions.

+-----------------+-----------------------------+
+ Name:           | ContourWriter               |
+-----------------+-----------------------------+
| Source:         | pyfem/io/ContourWriter.py   |
+-----------------+-----------------------------+

---------------------
Mandatory parameters:
---------------------

+-----------------+-----------------------------+
| Parameter name  | Description                 |
+-----------------+-----------------------------+
| ``nodes``       | An array of nodes.          |
+-----------------+-----------------------------+

------------
Examples:
------------

ch13: PeelTest40.pro

