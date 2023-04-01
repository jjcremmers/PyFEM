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

import logging

def setLogger( props ):
  
  level = "normal"
  
  if hasattr(props,"logger"):
    level = props.logger.level
    
    if level not in ["normal","info","debug","critical","warning","silent"]:
      raise NotImplementedError('Logger level should be "normal", "info", "debug", "critical", "silent" or "warning"')
    
  logger    = logging.getLogger()
  handler   = logging.StreamHandler()
    
  if level == "debug":
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    logger .setLevel(logging.DEBUG)
  elif level == "critical" or level == "silent":
    formatter = logging.Formatter('  %(message)s')
    logger .setLevel(logging.CRITICAL)    
  elif level == "warning":
    formatter = logging.Formatter('  %(message)s')
    logger .setLevel(logging.WARNING)      
  else:
    formatter = logging.Formatter('  %(message)s')
    logger .setLevel(logging.INFO)
    
  handler.setFormatter(formatter)
  logger .addHandler(handler)

  return logger

def getLogger():

  return logging.getLogger()
