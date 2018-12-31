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
from pyfem.util.dataStructures import Properties
from pyfem.util.dataStructures import GlobalData

from pyfem.fem.NodeSet    import NodeSet
from pyfem.fem.ElementSet import ElementSet
from pyfem.fem.DofSpace   import DofSpace

from pyfem.util.fileParser import fileParser

import getopt,os.path

def InputReader( argv ):

  options, remainder = getopt.getopt( argv , 'a:k:v', ['all','author='])

  proFileName  = argv[1]

  return InputRead( proFileName )

def InputRead( fname ):

  if fname[-4:] == '.pro':
    props        = fileParser( fname )
  else:
    props        = fileParser( fname+'.pro')

  dataFileName = props.input

  nodes = NodeSet()
  nodes.readFromFile( dataFileName )
  
  elems = ElementSet( nodes , props )
  elems.readFromFile( dataFileName )
  
  dofs = DofSpace( elems )
  dofs.readFromFile( dataFileName )

  globdat = GlobalData( nodes, elems, dofs ) 

  globdat.readFromFile( dataFileName )

  globdat.active = True
  globdat.cycle  = 0
  globdat.prefix = os.path.splitext(fname)[0]
	
  return props,globdat
