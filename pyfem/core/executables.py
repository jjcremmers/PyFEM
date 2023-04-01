################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by                #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since        #
#  then augmented and maintained by Joris J.C. Remmers.                        #
#  All rights reserved.                                                        #
#                                                                              #
#  A github repository, with the most up to date version of the code,          #
#  can be found here:                                                          #
#     https://github.com/jjcremmers/PyFEM/                                     #
#     https://pyfem.readthedocs.io/                                            #	
#                                                                              #
#  The original code can be downloaded from the web-site:                      #
#     http://www.wiley.com/go/deborst                                          #
#                                                                              #
#  The code is open source and intended for educational and scientific         #
#  purposes only. If you use PyFEM in your research, the developers would      #
#  be grateful if you could cite the book.                                     #    
#                                                                              #
#  Disclaimer:                                                                 #
#  The authors reserve all rights but do not guarantee that the code is        #
#  free from errors. Furthermore, the authors shall not be liable in any       #
#  event caused by the use of the program.                                     #
################################################################################

import sys,os
sys.path.insert(0, os.getcwd() )

from pyfem.io.InputReader   import InputRead
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver   import Solver

def run( fileName ):

  props,globdat = InputRead( fileName )

  solver = Solver        ( props , globdat )
  output = OutputManager ( props , globdat )

  while globdat.active:
    solver.run( props , globdat )
    output.run( props , globdat )

  print("PyFem analysis terminated successfully.")

def readInput( fileName ):

  props,globdat = InputRead( fileName )

  solver = Solver        ( props , globdat )
  output = OutputManager ( props , globdat )

  globdat.props = props

  return globdat

def calcSingleStep( globdat ):

  solver.run( globdat.props , globdat )
