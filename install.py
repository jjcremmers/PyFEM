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
import subprocess
import sys

import matplotlib
import numpy
import scipy

print("\n ===============================================================\n")

# get operating system

osName = sys.platform
osFullName = osName

if osName == "linux":
  osFullName = "Linux"
elif osName == "win32":
  osFullName = "Windows"

print(f"  Operating system                       : {osFullName}")

# check python version

versionLong = sys.version.split(' ')
version = versionLong[0].split('.')

print(f"  Python version detected     {versionLong[0]:>10s} : ", end=' ')

if int(version[0]) == 3 and int(version[1]) >= 6:
  print("   OK")
elif int(version[0]) == 2:
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

  versionLong = numpy.__version__
  version = versionLong.split('.')

  print(f"  Numpy version detected      {versionLong:>10s} : ", end=' ')

  if int(version[0]) == 1 and int(version[1]) >= 6:
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

  versionLong = scipy.__version__
  version = versionLong.split('.')

  print(f"  Scipy version detected      {versionLong:>10s} : ", end=' ')

  if int(version[0]) == 0 and int(version[1]) >= 9:
    print("   OK")
  elif int(version[0]) >= 1 and int(version[1]) >= 0:
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

  versionLong = matplotlib.__version__
  version = versionLong.split('.')

  print(f"  Matplotlib version detected {versionLong:>10s} : ", end=' ')

  if int(version[0]) >= 1 and int(version[1]) >= 0:
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

  versionLong = meshio.__version__
  version = versionLong.split('.')

  print(f"  Meshio version detected     {versionLong:>10s} : ", end=' ')

  if int(version[0]) <= 3:
    print("  Not OK\n")
    answer = input("    Do you want to install the latest version meshio? (Y/N)\n")
    if answer.lower() == "y" or answer.lower() == 'yes':
      subprocess.run(['pip', 'install', 'meshio'])
    else:
      print("\n    You cannot use gmsh input files!\n")
  else:
    print("   OK")
except ImportError:
  print("  Meshio not detected                    : Not OK")
  answer = input("    Do you want to install the latest version meshio? (Y/N)\n")
  if answer.lower() == "y" or answer.lower() == 'yes':
    subprocess.run(['pip', 'install', 'meshio'])
  else:
    print("\n    You cannot use gmsh input files!\n")

# check pickle version

try:
  import pickle

  versionLong = pickle.format_version
  version = versionLong.split('.')

  print(f"  Pickle version detected     {versionLong:>10s} : ", end=' ')

  if int(version[0]) >= 4:
    print("   OK")
except ImportError:
  print("  pickle not detected                    : Not OK")
  print("\n    Please install pickle\n")
  print("      'pip install pickle'\n")
  print("    or run PyFEM with limited functionality.\n")

# check h5py version

try:
  import h5py

  versionLong = h5py.__version__
  version = versionLong.split('.')

  print(f"  H5py version detected       {versionLong:>10s} : ", end=' ')

  if int(version[0]) >= 2:
    print("   OK")
except ImportError:
  print("  h5py not detected                    : Not OK")
  print("\n    Please install h5py\n")
  print("      'pip install h5py'\n")
  print("    or run PyFEM with limited functionality.\n")

# check PySide version

try:
  import PySide6

  versionLong = PySide6.__version__
  version = versionLong.split('.')

  print(f"  PySide version detected     {versionLong:>10s} : ", end=' ')

  if int(version[0]) >= 6 and int(version[1]) >= 0:
    print("   OK")
  else:
    print("    Please install PySide 6.0.0 or higher and reconfigure PyFEM.\n")
    sys.exit()
except ImportError:
  print("  PySide not detected                    : Not OK")
  print("\n    Please install PySide\n")
  print("      'pip install PySide'\n")
  print("    or run PyFEM with limited functionality.\n")

# check vtk version

try:
  import vtk

  versionLong = vtk.__version__
  version = versionLong.split('.')

  print(f"  vtk version detected        {versionLong:>10s} : ", end=' ')

  if int(version[0]) >= 9 and int(version[1]) >= 0:
    print("   OK")
  else:
    print("    Please install vtk 9.0.0 or higher and reconfigure PyFEM.\n")
    sys.exit()
except ImportError:
  print("  vtk not detected                       : Not OK")
  print("\n    Please install vtk 9.0.0 or higher\n")
  sys.exit()

# get current path

path = os.getcwd()

print("\n ===============================================================")
print("  INSTALLATION SUCCESSFUL!")
print(" ===============================================================\n")

if osName[:5] == "linux":

  batfile = open('pyfem.sh', 'w')

  fexec = sys.executable

  batfile.write(fexec + ' ' + path + '/PyFEM.py "$1"')

  batfile.close()

  subprocess.run(['chmod', '+x', 'pyfem.sh'])

  batfile = open('pyfem_gui.x', 'w')

  fexec = sys.executable

  batfile.write(fexec + ' ' + path + '/pyfem_gui.py &')

  batfile.close()

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

elif osName[:6] == "darwin":

  print(" Add the following line to ~/.bashrc :\n")
  # print('   export PYTHONPATH="' + path + '"')
  print("    alias pyfem='python3 " + path + "/PyFEM.py'\n")
  print(" ===============================================================\n")
  print("  Installation successful!")
  print("  See the user manual for further instructions.\n\n")

elif osName[:3] == "win":

  fexec = sys.executable

  if fexec[-5:] == "w.exe":
    fexec = fexec[:-5] + ".exe"

  batfile = open('pyfem.bat', 'w')
  batfile.write(fexec + ' ' + path + '\\PyFEM.py %1')
  batfile.close()

  batfile = open('pyfem.exe', 'w')
  batfile.write(fexec + ' ' + path + '\\pyfem_gui.py')
  batfile.close()

  print("  You can run PyFEM from any directory by typing:\n")
  print("    [path_to_this_directory]\\pyfem inputFile.pro\n")
  print("  or, when you add this directory to your system path,")
  print("  you can run it by typing:\n")
  print("    pyfem inputFile.pro\n")
  print("  See the user manual for further instructions.\n")

else:
  print("Operating system ", osName, " not known.")

sys.exit()
