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

from typing import Any, List
from pyfem.models.BaseModel import BaseModel


class ModelManager:
    """Manager for multiple model instances in PyFEM simulations.

    The ModelManager dynamically loads and coordinates multiple models based on
    the input file configuration. It handles model instantiation, property
    assignment, and orchestrates the execution of all active models during
    the simulation.

    Attributes
    ----------
    models : List[BaseModel]
        List of instantiated model objects (e.g., Contact, RVE).

    Notes
    -----
    Currently supported model types:
    - Contact: Contact mechanics models
    - RVE: Representative Volume Element models with periodic BCs
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """Initialize the ModelManager and instantiate all configured models.

        Parameters
        ----------
        props : Properties
            Global properties object. If a property exists with the same name as
            a supported model type, that model will be instantiated with those
            properties.
        globdat : GlobalData
            Global data structure containing nodes, elements, dofs, and state vectors.

        Notes
        -----
        The initialization dynamically imports and creates model instances based on
        the presence of corresponding attributes in props. For each model type found,
        the model class is imported from pyfem.models.<ModelType> and instantiated.
        """

        modelTypes = ["Contact", "RVE"]
        
        self.models: List[BaseModel] = []

        for modelType in modelTypes:
            if hasattr(props, modelType):
                modelProps = getattr(props, modelType)
                modelProps.modelType = modelType
                
                exec("from pyfem.models." + modelType + " import " + modelType)
                
                self.models.append(eval(modelType + "(modelProps, globdat)"))

    def run(self, props: Any, globdat: Any) -> None:
        """Execute all active models in sequence.

        Parameters
        ----------
        props : Properties
            Global properties object passed to each model's run method.
        globdat : GlobalData
            Global data structure containing nodes, elements, dofs, and state vectors.
            Models may modify constraints, forces, or other global data.

        Notes
        -----
        Models are executed in the order they were added during initialization.
        Each model's run() method is called, allowing it to apply constraints,
        forces, or other modifications to the global data structure.
        """

        for model in self.models:
            model.run(props, globdat)      
