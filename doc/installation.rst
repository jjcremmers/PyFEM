Installation
============

Below, instructions are given how to access PyFEM and install it on windows using Git Bash. 

The original code can be downloaded from the website that accompanies the book.

  http://www.wiley.com/go/deborst
  
The latest development version is available on Github:

  https://github.com/jjcremmers/PyFEM

This version of the PyFEM is written to work properly in combination with 
Python version 3.x. In addition, the code uses the modules numpy, scipy and
matplotlib. Installation guidelines are given for various operating systems.

Linux
-----

The ''python'' program and the modules ''numpy'', ''scipy'' and ''matplotlib''
are included in most common distributions of Linux and can be installed without any problems. In many
cases, different versions of ''python'' are offered. Please make sure that ''python'' version 3.6 or higher is
installed. In addition, the modules ''meshio'', ''pickle'' and ''h5py'' can be installed for additional functionality.

Execute the file ''install.py'' in the root directory ''pyfem''. In a terminal, one can type:

  python3 install.py

This script will check if the correct versions of Python and the various modules are available. In addition,
the total path to the executable is given. For your own convenience, you can add this to your ''.bashrc'' file:

  alias  pyfem="python <pyfemdir>/PyFEM.py"

When using csh or tcsh add the following line to ''.cshrc'' or ''.tcshrc'':

  alias  pyfem "python <pyfemdir>/PyFEM.py"

The main program ''pyfem'' can be run from the command prompt. For example, in order to run the
file ''StressWave20x20.pro'' in the directory ''examples/ch05'', simply type:

  pyfem StressWave20x20.pro
