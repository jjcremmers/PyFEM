################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2024. The code is written in 2011-2012 by                #
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

import time
from pyfem.util.logger    import getLogger
from pyfem.util.plotUtils import plotTime

logger = getLogger()


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


class BaseModule:

  def __init__ ( self, props ):

    self.isSolver = False
     
    if hasattr(props,'currentModule') and hasattr(props,props.currentModule):
      currentModule = props.currentModule
      
      if 'solver' in currentModule:
        self.isSolver = True
    elif hasattr(props,'currentModule'):
      currentModule = props.currentModule
      
      if 'solver' in currentModule:
        self.isSolver = True
              
    else:
      currentModule = self.__class__.__name__
      
      print(currentModule)
      if currentModule.endswith("olver") == "olver":
        currentModule = "solver"
        
        self.isSolver = True
    
    c = currentModule.split('.')  
    
    if len(c) == 1:
      if hasattr(props,currentModule):
        self.myProps = getattr(props,currentModule)

        for name,val in self.myProps:
          setattr( self, name, val )
    elif len(c) == 2:
      if hasattr(props,c[0]):
        p2 = getattr(props,c[0])
        if hasattr(p2,c[1]):
          self.myProps = getattr(p2,c[1])
          
          for name,val in self.myProps:
            setattr( self, name, val ) 
    else:
      print("NNOO")
          
    self.type = self.__class__.__name__
 
    '''
def get_nested_attr(obj, attr_path, default=None):
    """
    Access a nested attribute in an object hierarchy.

    Parameters:
    - obj: The base object
    - attr_path: A string of dot-separated attribute names (e.g., "bar.foo.y")
    - default: The value to return if any attribute in the chain doesn't exist

    Returns:
    - The value of the nested attribute, or the default if any attribute is missing
    """
    try:
        for attr in attr_path.split('.'):
            obj = getattr(obj, attr)
        return obj
    except AttributeError:
        return default

# Example usage:
# Assuming x.bar.foo.y exists
value = get_nested_attr(x, "bar.foo.y", default="Attribute not found")
print(value)
    '''

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    
  def writeHeader( self , cycle = None ):
  
    self.t0     = time.time()
    cycleString = ""
    
    if cycle != None:
      cycleString = "  step: " + str(cycle);
    
    if self.isSolver:
      logger.info("")
      logger.info("=============================================================")
      logger.info("  " + self.type + cycleString )
      logger.info("=============================================================")      
    else:
      logger.debug("-------------------------------------------------------------")    
      logger.debug("  Module " + self.type)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

      
  def writeFooter( self , globdat ):
  
    t1 = time.time()
    
    if self.isSolver:
      logger.info("    Elapsed time (this step).. : " + plotTime(t1-self.t0))
      logger.info("    Total elapsed time........ : " + plotTime(t1-globdat.startTime))    
    else:
      logger.debug("    Elapsed time (this step).. : " + plotTime(t1-self.t0))
      logger.debug("    Total elapsed time........ : " + plotTime(t1-globdat.startTime))         

