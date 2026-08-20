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

import os.path,pickle,time

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
    
  if parameters is not None:  
    for p in parameters:
      x = p.split("=")
      props.store(x[0],x[1])
      
  if dname is not None:
    return props,data["globdat"]

  pathName, _ = os.path.split(fname)

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
  from pyfem.core.cli import parse_arguments

  args = parse_arguments(argv)
  return args.input_file, args.dump_file, args.parameters
