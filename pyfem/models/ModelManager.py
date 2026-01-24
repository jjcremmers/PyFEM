# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers


from typing import Any, List
from importlib import import_module
from pyfem.models.BaseModel import BaseModel


class ModelManager():
    """
    Manager for multiple model instances in PyFEM simulations.

    Dynamically loads and coordinates multiple models based on the input file configuration.
    Handles model instantiation, property assignment, and orchestrates the execution
    of all active models during the simulation.

    Attributes:
        models (List[BaseModel]): List of instantiated model objects (e.g., Contact, RVE).
    """

    def __init__(self, props: Any, globdat: Any) -> None:
        """
        Initialize the ModelManager and instantiate all configured models.

        Args:
            props (Any): Global properties object. If a property exists with the same name as
                a supported model type, that model will be instantiated with those properties.
            globdat (Any): Global data structure containing nodes, elements, dofs, and state vectors.
        """
        modelTypes = ["Contact", "RVE"]
        self.models: List[BaseModel] = []

        for modelType in modelTypes:
            if hasattr(props, modelType):
                modelProps = getattr(props, modelType)
                modelProps.modelType = modelType
                try:
                    module = import_module(f"pyfem.models.{modelType}")
                except ModuleNotFoundError as e:
                    raise ImportError(
                        f"Model module 'pyfem.models.{modelType}' not found. "
                        f"Check the 'type' in your input file."
                    ) from e
                try:
                    model_cls = getattr(module, modelType)
                except AttributeError as e:
                    raise ImportError(
                        f"Class '{modelType}' not found in module 'pyfem.models.{modelType}'. "
                        f"Ensure the class name matches the file name."
                    ) from e
                self.models.append(model_cls(modelProps, globdat))

    def takeAction(self, action, mbuilder, props: Any, globdat: Any) -> None:
        """
        Execute all active models in sequence.

        Args:
            action: Action to perform (currently ignored, for compatibility).
            mbuilder: ModelBuilder instance (currently ignored, for compatibility).
            props (Any): Global properties object passed to each model's run method.
            globdat (Any): Global data structure containing nodes, elements, dofs, and state vectors.
        """
        
        for model in self.models:
            model.run(props, globdat)
