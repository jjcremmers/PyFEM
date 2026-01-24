# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers


from typing import Any, List
from importlib import import_module
from pyfem.models.BaseModel import BaseModel

from pyfem.util.logger import getLogger, separator, logVariable

logger = getLogger()

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
        self.modelss = []

        if hasattr(props, "models"):
            self._getModels(props, globdat)
        
    def _getModels(self, props: Any, globdat: Any) -> None:
        """
        Private method to load and instantiate models from props.models.

        Args:
            props (Any): Global properties object containing 'models' attribute.
            globdat (Any): Global data/state object.
        """

        self.models = props.models
        for model in self.models:
            modelProps = getattr(props, model)
            modelType = modelProps.type
            try:
                mode = import_module(f"pyfem.models.{modelType}")
            except ModuleNotFoundError as e:
                raise ImportError(
                    f"Model module 'pyfem.models.{modelType}' not found. "
                    f"Check the 'type' in your input file."
                ) from e
            try:
                model_cls = getattr(mode, modelType)
            except AttributeError as e:
                raise ImportError(
                    f"Class '{modelType}' not found in module 'pyfem.materials.{modelType}'. "
                    f"Ensure the class name matches the file name."
                ) from e
            self.modelss.append(model_cls(modelProps, globdat))

        separator("=")
        logger.info("Reading models")

        for model in self.models:
            logger.info(f"  {model}")      

        separator("=")
            
    def takeAction(self, action: str, mbuilder, props: Any, globdat: Any) -> None:
        """
        Execute the method named by the string 'action' on all active models, if it exists.

        Args:
            action (str): Name of the method to call on each model (as a string).
            mbuilder: MatrixBuilder instance to pass as argument.
            props (Any): Global properties object passed to the model method.
            globdat (Any): Global data structure containing nodes, elements, dofs, and state vectors.
        """

        for model in self.modelss:
            method = getattr(model, action, None)
            if callable(method):
                method(props, globdat, mbuilder)
