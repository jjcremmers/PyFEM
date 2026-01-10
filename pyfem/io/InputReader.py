# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util.dataStructures import Properties
from pyfem.util.dataStructures import GlobalData

from pyfem.fem.NodeSet         import NodeSet
from pyfem.fem.ElementSet      import ElementSet
from pyfem.fem.DofSpace        import DofSpace

from pyfem.models.ModelManager import ModelManager

from pyfem.util.fileParser     import fileParser
from pyfem.util.logger         import setLogger, separator

import getopt,os.path,pickle,time

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

  t1 = time.time()
  
  if dname is not None:
    with open(dname, 'rb') as f:
      data = pickle.load(f)
      props = data["props"]
  
  if fname is not None:
    if fname[-4:] == '.pro':
      props        = fileParser( fname )
    else:
      props        = fileParser( fname+'.pro')
    
  pathName, _ = os.path.split(fname)
  
  if parameters is not None:  
    for p in parameters:
      x = p.split("=")
      props.store(x[0],x[1])
      
  if dname is not None:
    return props,data["globdat"]

  dataFileName = props.input
  
  dataFileName = os.path.join(pathName,dataFileName)
 
  logger = setLogger( props )
  
  separator("=")
  logger.info("  PyFEM analysis: " + fname )
  separator("=")

  nodes = NodeSet()
  nodes.readFromFile( dataFileName )
  
  elems = ElementSet( nodes , props )
  elems.readFromFile( dataFileName )
  elems.logInfo()
  
  dofs = DofSpace( elems )
  dofs.readFromFile( dataFileName )

  globdat = GlobalData( nodes, elems, dofs ) 

  globdat.readFromFile( dataFileName )

  globdat.active = True
  globdat.prefix = os.path.splitext(fname)[0]
   
  globdat.models  = ModelManager( props , globdat )
  
  globdat.startTime = t1  
  	
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
