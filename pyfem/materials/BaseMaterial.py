############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stable version can be downloaded from the web-site:          #
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/jjcremmers/PyFEM                                  #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #  
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################

import copy
from typing import Any, Dict, List, Tuple, Union
import numpy as np


class BaseMaterial:
    """
    Base class for material models in finite element analysis.
    
    This class provides the fundamental structure and methods for implementing
    constitutive material models. It handles material properties, history
    variables for path-dependent materials, and output data management.
    
    Attributes
    ----------
    numericalTangent : bool
        Flag indicating whether to use numerical tangent computation.
    storeOutputFlag : bool
        Flag indicating whether to store output data.
    oldHistory : Dict[str, Any]
        Dictionary storing history variables from the previous converged step.
    newHistory : Dict[str, Any]
        Dictionary storing history variables for the current step.
    outLabels : List[str]
        Labels for output variables.
    outData : ndarray
        Array containing output data values.
    solverStat : Any
        Solver statistics object.
    """

    def __init__(self, props) -> None:
        """
        Initialize the BaseMaterial instance.
        
        Parameters
        ----------
        props : object
            Properties object containing material parameters and solver statistics.
            The props object should be iterable as (name, value) pairs and have
            a solverStat attribute.
        """
        self.numericalTangent = False
        self.storeOutputFlag = False
        
        for name, val in props:
            setattr(self, name, val)

        self.oldHistory = {}
        self.newHistory = {}

        self.outLabels = []
        self.solverStat = props.solverStat
        
    def setHistoryParameter(self, name: str, val: Any) -> None:
        """
        Set a history parameter for the current step.
        
        History parameters are used to store internal state variables for
        path-dependent material models (e.g., plastic strain, damage variables).
        
        Parameters
        ----------
        name : str
            Name of the history parameter.
        val : Any
            Value of the history parameter (typically float or ndarray).
        """
        self.newHistory[name] = val
        return
           
    def getHistoryParameter(self, name: str) -> Union[float, np.ndarray]:
        """
        Retrieve a history parameter from the previous converged step.
        
        Parameters
        ----------
        name : str
            Name of the history parameter to retrieve.
        
        Returns
        -------
        Union[float, ndarray]
            The value of the history parameter. Returns a copy for array types
            to prevent unintended modifications.
        """
        if type(self.oldHistory[name]) == float:
            return self.oldHistory[name]
        else:
            return self.oldHistory[name].copy()
                  
    def commitHistory(self) -> None:
        """
        Commit the current history variables to old history.
        
        This method is called when a load step has converged, copying the
        current (new) history variables to become the old history variables
        for the next step. Uses deep copy to ensure complete independence.
        """
        self.oldHistory = copy.deepcopy(self.newHistory)
        
    def setOutputLabels(self, labels: List[str]) -> None:
        """
        Set the labels for output variables.
        
        Initializes the output data array with zeros based on the number
        of output labels.
        
        Parameters
        ----------
        labels : List[str]
            List of strings representing the names of output variables
            (e.g., ['stress_xx', 'stress_yy', 'plastic_strain']).
        """
        self.outLabels = labels
        self.outData = np.zeros(len(self.outLabels))
        return
        
    def storeOutputs(self, data: np.ndarray) -> None:
        """
        Store output data if the store output flag is enabled.
        
        Parameters
        ----------
        data : ndarray
            Array of output data values corresponding to the output labels.
        """
        if self.storeOutputFlag:
            self.outData = data
        return
  
