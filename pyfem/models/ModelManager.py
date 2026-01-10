# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

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
