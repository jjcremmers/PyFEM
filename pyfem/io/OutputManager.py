################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
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

from importlib import import_module

class OutputManager:

    def __init__( self , props , globdat ):

        self.outman = []

        outputModules = props.outputModules

        for name in outputModules:
   
            props.currentModule = name

            outputType = name

            if hasattr( props , name):
                moduleProps = getattr( props, name )
                if hasattr( moduleProps , "type" ):
                    outputType = moduleProps.type

                try:
                    io = import_module(f"pyfem.io.{outputType}")
                    outputClass = getattr(io, outputType)  
                except ModuleNotFoundError as e:
                    raise ImportError(
                        f"Solver module 'pyfem.io.{outputType}' not found."
                    ) from e
                except AttributeError as e:
                    raise ImportError(
                        f"Class '{outputType}' not found in module "
                        f"'pyfem.io.{outputType}'."
                    ) from e

            self.outman.append( outputClass(props, globdat) )

            #exec("from pyfem.io."+ioType+" import "+ioType)

            #self.outman.append(eval(ioType+"( props , globdat )"))

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def run( self , props , globdat ):

        for i,output in enumerate(self.outman):
            output.run( props , globdat )
