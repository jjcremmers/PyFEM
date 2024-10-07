############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  Copyright (C) 2011-2024. The code is written in 2011-2012 by            #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since    #
#  then augmented and  maintained by Joris J.C. Remmers.                   #
#  All rights reserved.                                                    #
#                                                                          #
#  The latest stable version can be downloaded from the web-site:          #
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/jjcremmers/PyFEM                                  #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################

import os
import pickle
import platform
import subprocess
import sys

def _parse_version_string(version_string):
  # [:3] is used to ignore any additional version information
  v = version_string.split('.')[:3]
  return '.'.join(v), tuple(map(int, v))

print("\n ===============================================================\n")

# get operating system

os_name = platform.system()
print(f"  Operating system            {os_name:>10s} : ", end=' ')
if os_name in ("Linux", "Darwin", "Windows"):
  print("   OK")
else:
  print("  Not OK\n\n")
  print("    PyFEM is not supported on this operating system.\n")
  sys.exit()

# check python version

version_long, version = _parse_version_string(platform.python_version())

print(f"  Python version detected     {version_long:>10s} : ", end=' ')

if version >= (3, 6):
  print("   OK")
elif version[0] == 2:
  print("  Please note that PyFEM has been migrated to Python 3.x\n")
  print("    Install the latest version of Python 3.x and reconfigure PyFEM.\n")
  sys.exit()
else:
  print("  Not OK\n\n")
  print("    Install the latest version of Python 3.x and reconfigure PyFEM.\n")
  sys.exit()

# check numpy version

try:
  import numpy

  version_long, version = _parse_version_string(numpy.__version__)

  print(f"  Numpy version detected      {version_long:>10s} : ", end=' ')

  if version >= (1, 6):
    print("   OK")
  else:
    print("  Not OK\n\n")
    print("    Please install Numpy 1.6.x or higher and reconfigure PyFEM.\n")
    sys.exit()
except ImportError:
  print("  NumPy not detected                      : Not OK")
  print("    Please install Numpy 1.6.x or higher and reconfigure PyFEM.\n")
  sys.exit()

# check scipy version

try:
  import scipy

  version_long, version = _parse_version_string(scipy.__version__)

  print(f"  Scipy version detected      {version_long:>10s} : ", end=' ')

  if version >= (0, 9):
    print("   OK")
  else:
    print("    Please install Scipy 0.9.x or higher and reconfigure PyFEM.\n")
    sys.exit()
except ImportError:
  print("  SciPy not detected                     : Not OK")
  print("    Please install Scipy 0.9.x or higher and reconfigure PyFEM.\n")
  sys.exit()

# check matplotlib

try:
  import matplotlib

  version_long, version = _parse_version_string(matplotlib.__version__)

  print(f"  Matplotlib version detected {version_long:>10s} : ", end=' ')

  if version >= (1, 0):
    print("   OK")
  else:
    print("  Not OK\n\n   Please install Matplotlib 1.0.x or higher\n")
except ImportError:
  print("  matplotlib not detected                : Not OK")
  print("\n    Please install Matplotlib 1.0.x or higher\n")
  sys.exit()

# check meshio version

try:
  import meshio

  version_long, version = _parse_version_string(meshio.__version__)

  print(f"  Meshio version detected     {version_long:>10s} : ", end=' ')

  if version >= (4, 0):
    print("   OK")
  else:
    print("  Not OK\n")
    answer = input("    Do you want to install the latest version meshio? (Y/N)\n")
    if answer.lower() == "y" or answer.lower() == 'yes':
      subprocess.run(['pip', 'install', 'meshio'], check=True)
    else:
      print("\n    You cannot use gmsh input files!\n")
except ImportError:
  print("  Meshio not detected                    : Not OK")
  answer = input("    Do you want to install the latest version meshio? (Y/N)\n")
  if answer.lower() == "y" or answer.lower() == 'yes':
    subprocess.run(['pip', 'install', 'meshio'], check=True)
  else:
    print("\n    You cannot use gmsh input files!\n")

# check h5py version

try:
  import h5py

  version_long, version = _parse_version_string(h5py.__version__)

  print(f"  H5py version detected       {version_long:>10s} : ", end=' ')

  if version >= (2, 0):
    print("   OK")
except ImportError:
  print("  h5py not detected                      : Not OK")
  answer = input("    Do you want to install the latest version of h5py? (Y/N)\n")
  if answer.lower() == "y" or answer.lower() == 'yes':
    subprocess.run(['pip', 'install', 'h5py'], check=True)
  else:
    print("\n    You cannot write h5 files!\n")

# check PySide version

try:
  import PySide6

  version_long, version = _parse_version_string(PySide6.__version__)

  print(f"  PySide version detected     {version_long:>10s} : ", end=' ')

  if version >= (6, 0):
    print("   OK")
  else:
    print("    Please install PySide 6.0.0 or higher and reconfigure PyFEM.\n")
    sys.exit()
except ImportError:
  print("  PySide not detected                    : Not OK")
  answer = input("    Do you want to install the latest version of PySide6? (Y/N)\n")
  if answer.lower() == "y" or answer.lower() == 'yes':
    subprocess.run(['pip', 'install', 'PySide6'], check=True)
  else:
    print("    or run PyFEM with limited functionality.\n")

# check vtk version

try:
  import vtk

  version_long, version = _parse_version_string(vtk.__version__)

  print(f"  vtk version detected        {version_long:>10s} : ", end=' ')

  if version >= (9, 0):
    print("   OK")
  else:
    print("    Please install vtk 9.0.0 or higher and reconfigure PyFEM.\n")
    sys.exit()
except ImportError:
  print("  vtk not detected                       : Not OK")
  answer = input("    Do you want to install the latest version of vtk? (Y/N)\n")
  if answer.lower() == "y" or answer.lower() == 'yes':
    subprocess.run(['pip', 'install', 'vtk'], check=True)
  else:
    print("    or run PyFEM with limited functionality.\n")

# get current path

path = os.getcwd()

print("\n ===============================================================")
print("  INSTALLATION SUCCESSFUL!")
print(" ===============================================================\n")

if os_name == "Linux":

  with open('pyfem.sh', 'w') as bat_file:
    fexec = sys.executable
    bat_file.write(fexec + ' ' + path + '/PyFEM.py "$1"')

  subprocess.run(['chmod', '+x', 'pyfem.sh'])

  with open('pyfem_gui.x', 'w') as bat_file:
    fexec = sys.executable
    bat_file.write(fexec + ' ' + path + '/pyfem_gui.py &')

  subprocess.run(['chmod', '+x', 'pyfem_gui.x'])

  print("  You can run PyFEM in command line from any directory by typing:\n")
  print("    [relative_path_to_this_directory]/pyfem.sh inputFile.pro\n")
  print("  The gui can be run by typing from any directory:\n")
  print("    [relative_path_to_this_directory]/pyfem_gui.exe\n")
  print("  Alternatively, you can make an aliases. When using a bash shell,")
  print("  add the following lines to the file ~/.bashrc :\n")
  print("   alias pyfem='python3 " + path + "/PyFEM.py'")
  print("   alias pyfem_gui = '" + path + "/pyfem_gui.x'\n")
  print("  and you can run PyFEM in commandline from any directory by typing:\n")
  print("    pyfem inputFile.pro\n")
  print("  and the gui by typing:\n")
  print("    pyfem_gui\n")
  print("  See the user manual for further instructions.\n")

elif os_name == "Darwin":

  print(" Add the following line to ~/.bashrc :\n")
  # print('   export PYTHONPATH="' + path + '"')
  print("    alias pyfem='python3 " + path + "/PyFEM.py'\n")
  print(" ===============================================================\n")
  print("  Installation successful!")
  print("  See the user manual for further instructions.\n\n")

elif os_name == "Windows":

  fexec = sys.executable

  if fexec[-5:] == "w.exe":
    fexec = fexec[:-5] + ".exe"
    
  with open('pyfem.bat', 'w') as bat_file:
    bat_file.write(fexec + ' ' + path + '\\PyFEM.py %1')

  with open('pyfem_gui.exe', 'w') as bat_file:
    bat_file.write(fexec + ' ' + path + '\\pyfem_gui.py')    

  print("  You can run PyFEM from any directory by typing:\n")
  print("    [path_to_this_directory]\\pyfem inputFile.pro\n")
  print("  or, when you add this directory to your system path,")
  print("  you can run it by typing:\n")
  print("    pyfem inputFile.pro\n")
  print("  See the user manual for further instructions.\n")

sys.exit()
