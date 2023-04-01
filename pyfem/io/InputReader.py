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

from pyfem.util.dataStructures import Properties
from pyfem.util.dataStructures import GlobalData

from pyfem.fem.NodeSet     import NodeSet
from pyfem.fem.ElementSet  import ElementSet
from pyfem.fem.DofSpace    import DofSpace
from pyfem.fem.Contact     import Contact

from pyfem.util.fileParser import fileParser
from pyfem.util.logger     import setLogger

import getopt,os.path,pickle

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def InputReader( argv ):

  pName,dName,params = getArguments( argv )
  
  return InputRead( pName , dName , params )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def InputRead( fname , dname = None , parameters = None ):

  if dname is not None:
    with open(dname, 'rb') as f:
      data = pickle.load(f)
      props = data["props"]
  
  if fname is not None:
    if fname[-4:] == '.pro':
      props        = fileParser( fname )
    else:
      props        = fileParser( fname+'.pro')
    
  if parameters is not None:  
    for p in parameters:
      x = p.split("=")
      props.store(x[0],x[1])
      
  if dname is not None:
    return props,data["globdat"]

  dataFileName = props.input

  logger = setLogger( props )

  nodes = NodeSet()
  nodes.readFromFile( dataFileName )
  logger.info(nodes)
  
  elems = ElementSet( nodes , props )
  elems.readFromFile( dataFileName )
  logger.info(elems)  
  
  dofs = DofSpace( elems )
  dofs.readFromFile( dataFileName )

  globdat = GlobalData( nodes, elems, dofs ) 

  globdat.readFromFile( dataFileName )

  globdat.active = True
  globdat.prefix = os.path.splitext(fname)[0]
  
  globdat.contact = Contact( props )
  	
  return props,globdat
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def getArguments( argv ):

  slist = 'd:i:hvp:'
  llist = ['dump=','input=','help','version']
  
  options, remainder = getopt.getopt( argv[1:] , slist , llist )

  proFileName  = None
  dumpFileName = None
  parameters   = []
  
  if len(options) == 0:
    proFileName  = argv[1]
    options, remainder = getopt.getopt( argv[2:] , slist, llist )
    
  for opt, arg in options:      
    if opt in ('-i', '--input'):
      proFileName  = arg
    elif opt in ('-d', '--dump'):
      dumpFileName  = arg
    elif opt in ('-h', '--help'):
      print("Help")
    elif opt in ('-p' , '--param'):
      parameters.append(arg)
      
  return proFileName,dumpFileName,parameters      
