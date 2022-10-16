********************************
Explicit time integration solver
********************************

The explicit time integration solver is discussed in detail in 
Section 5.2 of the book. The source code is explained in detail in Section 
5.3.

===========
Input
===========

+--------+--------------------------------+
|Name:   | ExplicitSolver                 |
+--------+--------------------------------+
|Source: | pyfem/solver/ExplicitSolver.py |
+--------+--------------------------------+

Mandatory parameters:

+--------+--------------------------------+
|dtime   | Magnitude of time step         |
+--------+--------------------------------+
|lam     | Load factor λ as a function of |
|        | time. This can be written as a |
|        | string. For example,           |
|        | ’4.0*sin(3.0*t)’ represents    |
|        | a sinusoidal load, with period |
|        | 3.0 and amplitude 4.0.         |
+--------+--------------------------------+

Optional parameters:

+---------+--------------------------------+
|maxCycle | Number of cycles after which   |
|         | the simulation will be         |
|         | terminated.                    |
+---------+--------------------------------+
|maxTime  | Time after which the           |
|         | simulation will be terminated. |
+---------+--------------------------------+

Examples:

ch05/StressWave20x20.pro
