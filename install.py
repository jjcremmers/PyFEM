############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stabke version can be downloaded from the web-site:          #
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
import os,sys,numpy,scipy,matplotlib

from PyQt5.Qt import PYQT_VERSION_STR

print("\n ===============================================================\n")

# get operating system

osName = sys.platform

# check python version

versionLong = sys.version.split(' ')
version     = versionLong[0].split('.')

print(" Python version detected     %10s : " %(versionLong[0]) , end=' '  )

if int(version[0]) == 3 and int(version[1]) >= 6:
  print(" OK")
elif int(version[0]) == 2:
  print(" Please note that PyFEM has been migrated to Python 3.x\n")
  print("   Install Pyhon 3.x\n")
else:
  print(" Not OK\n\n   Please install Python 2.6.x or 2.7.x\n")
  
# check numpy version

versionLong = numpy.__version__
version     = versionLong.split('.')

print(" Numpy version detected      %10s : " %(versionLong) , end=' '  )

if int(version[0]) == 1 and int(version[1]) >= 6:
  print(" OK")
else:
  print(" Not OK\n\n   Please install Numpy 1.6.x or higher\n")

# check scipy version

versionLong = scipy.__version__
version     = versionLong.split('.')

print(" Scipy version detected      %10s : " %(versionLong) , end=' '  )

if int(version[0]) == 0 and int(version[1]) >= 9:
  print(" OK")
elif int(version[0]) >= 1 and int(version[1]) >= 0:
  print(" OK")
else:
  print(" Not OK\n\n   Please install Scipy 0.9.x or higher\n")

versionLong = matplotlib.__version__
version     = versionLong.split('.')

print(" Matplotlib version detected %10s : " %(versionLong) , end=' '  )

if int(version[0]) >= 1 and int(version[1]) >= 0:
  print(" OK")
else:
  print(" Not OK\n\n   Please install Matplotlib 1.0.x or higher\n")
  
versionLong = PYQT_VERSION_STR
version     = versionLong.split('.')

print(" PyQt version detected       %10s : " %(versionLong) , end=' '  )

if int(version[0]) >= 5:
  print(" OK")
else:
  print(" Not OK\n\n   Please install PyQt 5.x or higher\n")

# get current path

path = os.getcwd()

if osName[:5] == "linux":

  print("\n LINUX INSTALLATION")
  print(" ===============================================================\n")
  print(" When using a bash shell, add the following lines")
  print(" to ~/.bashrc :\n")
  print('   export PYTHONPATH="'+path+'"')
  print("   alias  pyfem='python3 "+path+"/PyFEM.py'\n")
  print(" When using csh or tcsh add the following lines to")
  print(" ~/.cshrc or ~/.tcshrc :\n")
  print("   setenv PYTHONPATH "+path)
  print("   alias  pyfem 'python3 "+path+"/PyFEM.py'\n")

elif osName[:6] == "darwin":

  print("\n MAC-OS INSTALLATION")
  print(" ===============================================================\n")
  print(" Add the following lines to ~/.bashrc :\n")
  print('   export PYTHONPATH="'+path+'"')
  print("   alias  pyfem='python3 "+path+"/PyFEM.py'\n")

elif osName[:3] == "win":

  batfile = open( 'pyfem.bat' , 'w' )

  fexec = sys.executable
  
  if fexec[-5:] == "w.exe":
    fexec = fexec[:-5] + ".exe"
    
  print(fexec)
  batfile.write(fexec+' '+path+'\PyFEM.py %1')

  batfile.close()

  print("\n WINDOWS INSTALLATION")
  print(" ===============================================================\n")
  print(" Add the following path to PYTHONPATH and PATH:\n")
  print("   ",path,"\n")
  print(" See the user manual for instructions\n\n")

else:
  print("Operating system ",osName," not known.")

input("Press Enter to continue...")

