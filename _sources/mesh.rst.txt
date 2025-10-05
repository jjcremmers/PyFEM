###################
Mesh output writer
###################

The mesh output writer saves all data during a simulation to the disk. 
The data is organised as follows: during a simulation, a single output 
file filename.pvd will be created which refers to the output of single 
steps, which are stored in the file filename-xx.vtu, where xx indicates 
the step number. This data can be visualised by opening the file 
filename.pvd in the external program Paraview.

+---------------+-------------------------+
| Name:         |   MeshWriter            |
+---------------+-------------------------+
| Source:       | pyfem/io/MeshWriter.py  |
+---------------+-------------------------+

----------------------
Mandatory parameters:
----------------------

None. By default, every step is stored in the files with
the prefix of the current job.

-------------------------
Optional parameters:
-------------------------

+---------------+--------------------------+
| prefix        | The prefix of the output |
|               | filename that will be    |
|               | used. By default, the    |
|               | prefix of the input      |
|               | filename is used.        |
+---------------+--------------------------+
| interval      | The interval (number of  |
|               | cycles) for which output |
|               | is stored. By default,   |
|               | every step is stored.    |
+---------------+--------------------------+
| elementgroup  | When specified, only the |
|               | elements in this group   |
|               | will be stored. By       |
|               | default, all elements    |
|               | will be stored.          |
+---------------+--------------------------+

Examples:
ch03: cantilever8.pro
ch05: StressWave20x20.pro
ch06: ContDamExample.pro
ch13: PeelTest.pro

