# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

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
