# Developing Element Classes

This guide explains how to implement new finite element formulations in PyFEM.
Element classes define the computational core of finite element analysis,
including shape functions, strain-displacement relations, and element-level
matrix assembly.

## Overview

Elements in PyFEM inherit from the `Element` base class and implement methods
for computing:

- Tangent stiffness matrices
- Internal force vectors
- Mass matrices (for dynamic analysis)
- Output quantities at integration points

Each element type is responsible for:

1. Defining its degrees of freedom (`dofTypes`)
2. Computing strain-displacement matrices (B-matrices)
3. Calling material models to obtain stresses and tangent moduli
4. Integrating element matrices using numerical quadrature
5. Providing output data for post-processing

## Element Class Structure

All element classes must inherit from the `Element` base class located in
`pyfem/elements/Element.py`. The base class provides common functionality
for managing element data, material models, and output.

### Required Methods

An element class must implement at minimum:

```python
class MyElement(Element):

    def __init__(self, elnodes, props):
        """Initialize element with node IDs and properties."""
        Element.__init__(self, elnodes, props)
        # Define DOF types and dimensions
        self.dofTypes = ['u', 'v']  # For 2D problems
        
    def getTangentStiffness(self, elemdat):
        """Compute element stiffness matrix and internal forces."""
        # Implementation here
        
    def getInternalForce(self, elemdat):
        """Compute internal force vector only."""
        # Implementation here
```

### Optional Methods

For specialized functionality:

```python
def getMassMatrix(self, elemdat):
    """Compute element mass matrix for dynamic analysis."""
    pass
    
def getDissipation(self, elemdat):
    """Compute dissipated energy for arc-length methods."""
    pass
```

## Implementation Example

Here's a complete example implementing a 2D small-strain continuum element,
following the formulation presented in Chapter 4 of the book *"Non-Linear 
Finite Element Analysis of Solids and Structures"* by de Borst et al.

### Basic Structure

```python
# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Your Name

from .Element import Element
from pyfem.util.shapeFunctions import getElemShapeData
from pyfem.util.kinematics import Kinematics
import numpy as np

class MySmallStrainElement(Element):

    def __init__(self, elnodes, props):
        """Initialize the element.
        
        Args:
            elnodes: List of node IDs for this element
            props: Properties dictionary containing element parameters
        """
        Element.__init__(self, elnodes, props)
        
        # Define problem dimension from properties
        self.rank = props.rank  # 2 for 2D, 3 for 3D
        
        # Define degree of freedom types
        if self.rank == 2:
            self.dofTypes = ['u', 'v']
            self.nstr = 3  # Number of stress/strain components
        elif self.rank == 3:
            self.dofTypes = ['u', 'v', 'w']
            self.nstr = 6
            
        # Initialize kinematics utility
        self.kin = Kinematics(self.rank, self.nstr)
```

### Tangent Stiffness Method

The `getTangentStiffness` method computes the element stiffness matrix and
internal force vector. Following equation (4.42) in the book:

$$
\mathbf{K}^e = \int_{\Omega^e} \mathbf{B}^T \mathbf{D} \mathbf{B} \, d\Omega
$$

$$
\mathbf{f}_{int}^e = \int_{\Omega^e} \mathbf{B}^T \boldsymbol{\sigma} \, d\Omega
$$

```python
def getTangentStiffness(self, elemdat):
    """Compute element stiffness and internal forces.
    
    Args:
        elemdat: Element data container with:
            - coords: Nodal coordinates
            - state: Current displacement vector
            - Dstate: Displacement increment
            - stiff: Element stiffness matrix (output)
            - fint: Internal force vector (output)
    """
    # Get shape function data and integration points
    sData = getElemShapeData(elemdat.coords)
    
    # Loop over integration points
    for iData in sData:
        # Compute strain-displacement matrix (B-matrix)
        b = self.getBmatrix(iData.dhdx)
        
        # Compute strains from displacements
        self.kin.strain = b @ elemdat.state
        self.kin.dstrain = b @ elemdat.Dstate
        
        # Call material model to get stress and tangent
        sigma, tang = self.mat.getStress(self.kin)
        
        # Assemble element matrices (numerical integration)
        elemdat.stiff += b.T @ (tang @ b) * iData.weight
        elemdat.fint += b.T @ sigma * iData.weight
        
        # Store output data (stress, strain, etc.)
        self.appendNodalOutput(self.mat.outLabels(), 
                               self.mat.outData())
```

### Internal Force Method

For explicit dynamics or when only forces are needed:

```python
def getInternalForce(self, elemdat):
    """Compute internal force vector only (no stiffness).
    
    This method is more efficient when tangent stiffness is not needed,
    such as in explicit time integration schemes.
    """
    sData = getElemShapeData(elemdat.coords)
    
    for iData in sData:
        b = self.getBmatrix(iData.dhdx)
        
        self.kin.strain = b @ elemdat.state
        self.kin.dstrain = b @ elemdat.Dstate
        
        sigma, tang = self.mat.getStress(self.kin)
        elemdat.fint += b.T @ sigma * iData.weight
        
        self.appendNodalOutput(self.mat.outLabels(),
                               self.mat.outData())
```

### B-Matrix Computation

The strain-displacement matrix (B-matrix) relates nodal displacements to 
strains. For 2D small strains (from Section 4.3 of the book):

```python
def getBmatrix(self, dhdx):
    """Compute strain-displacement matrix.
    
    For 2D problems with strains [εxx, εyy, γxy]^T:
    
    Args:
        dhdx: Shape function derivatives [∂N/∂x, ∂N/∂y]
        
    Returns:
        b: Strain-displacement matrix (nstr × ndof)
    """
    b = np.zeros((self.nstr, len(self.dofTypes) * dhdx.shape[0]))
    
    if self.rank == 2:
        for i, dN in enumerate(dhdx):
            # εxx = ∂u/∂x
            b[0, i * 2] = dN[0]
            # εyy = ∂v/∂y  
            b[1, i * 2 + 1] = dN[1]
            # γxy = ∂u/∂y + ∂v/∂x
            b[2, i * 2] = dN[1]
            b[2, i * 2 + 1] = dN[0]
            
    elif self.rank == 3:
        # Similar for 3D case
        for i, dN in enumerate(dhdx):
            b[0, i * 3] = dN[0]     # εxx
            b[1, i * 3 + 1] = dN[1] # εyy
            b[2, i * 3 + 2] = dN[2] # εzz
            b[3, i * 3 + 1] = dN[2] # γyz
            b[3, i * 3 + 2] = dN[1]
            b[4, i * 3] = dN[2]     # γxz
            b[4, i * 3 + 2] = dN[0]
            b[5, i * 3] = dN[1]     # γxy
            b[5, i * 3 + 1] = dN[0]
            
    return b
```

### Mass Matrix (Optional)

For dynamic analysis, implement mass matrix computation:

```python
def getMassMatrix(self, elemdat):
    """Compute consistent element mass matrix.
    
    Following equation (9.1) in the book:
    M^e = ∫ ρ N^T N dΩ
    """
    sData = getElemShapeData(elemdat.coords)
    rho = self.rho  # Material density
    
    for iData in sData:
        N = self.getNmatrix(iData.h)
        elemdat.mass += rho * (N.T @ N) * iData.weight
        
def getNmatrix(self, h):
    """Shape function matrix for mass computation."""
    ndof = len(self.dofTypes)
    nnodes = len(h)
    N = np.zeros((ndof, ndof * nnodes))
    
    for i, hi in enumerate(h):
        for j in range(ndof):
            N[j, i * ndof + j] = hi
            
    return N
```

## Registration and Usage

### File Location

Place your new element class in the appropriate directory:

- `pyfem/elements/` for general elements
- Choose a descriptive filename, e.g., `MySmallStrainElement.py`

### Import in __init__.py

Add your element to `pyfem/elements/__init__.py`:

```python
from .MySmallStrainElement import MySmallStrainElement

__all__ = [
    'MySmallStrainElement',
    # ... other elements
]
```

### Using in Input Files

Users can now reference your element in `.pro` files:

```
MyElemGroup =
{
  type = "MySmallStrainElement";
  
  material = 
  {
    type = "PlaneStress";
    E    = 210.0e3;
    nu   = 0.3;
  };
};
```

## Testing and Validation

### Patch Tests

Always verify your element implementation with standard patch tests:

1. **Constant stress patch test** - Element should reproduce uniform stress
2. **Pure bending test** - For beam/plate elements
3. **Convergence studies** - Verify optimal convergence rates

Example test structure:

```python
# In test/testMyElement.py
import unittest
from pyfem.io.InputReader import InputRead

class TestMyElement(unittest.TestCase):

    def test_patch_test(self):
        """Verify element passes patch test."""
        props, globdat = InputRead("test/patch_test.pro")
        # Run analysis
        # Check if computed displacement matches theory
```

## Integration with Existing Code

### Material Interface

Your element should call material models through the standard interface:

```python
# Material returns stress and tangent modulus
sigma, tang = self.mat.getStress(self.kin)

# Material provides output labels and data
labels = self.mat.outLabels()
data = self.mat.outData()
```

### Shape Functions

Use the utility functions in `pyfem/util/shapeFunctions.py`:

```python
from pyfem.util.shapeFunctions import getElemShapeData

# Automatically selects quadrature based on element topology
sData = getElemShapeData(elemdat.coords)

# Returns integration point data with:
# - h: shape functions
# - dhdx: shape function derivatives
# - weight: integration weight
```

### Output Data

Provide element output for visualization:

```python
# Store nodal output (averaged to nodes)
self.appendNodalOutput(labels, data)

# Store element output (per element)
self.appendElementOutput(labels, data)
```

## Best Practices

1. **Follow naming conventions**: Use descriptive class names
2. **Document thoroughly**: Include docstrings for all methods
3. **Reference the book**: Cite equation numbers from de Borst et al.
4. **Test rigorously**: Implement patch tests and benchmark problems
5. **Handle edge cases**: Check for singular matrices, invalid input
6. **Optimize carefully**: Profile before optimizing
7. **Use utilities**: Leverage existing shape function and assembly routines

## Common Pitfalls

1. **Inconsistent coordinate systems**: Always use global coordinates
2. **Integration order**: Ensure sufficient quadrature points
3. **Locking phenomena**: Use selective integration or enhanced modes if needed
4. **Sign conventions**: Follow the book's sign conventions consistently
5. **Units**: Be consistent with units in material properties

## References

The theoretical foundation for element implementation can be found in:

*"Non-Linear Finite Element Analysis of Solids and Structures"*
by R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
John Wiley & Sons, 2012, ISBN 978-0470666449

Key chapters:

- Chapter 4: Isoparametric Elements
- Chapter 5: Element Technology
- Chapter 9: Dynamics and Time-Dependent Problems

## See Also

- [materials_dev.md](materials_dev.md) - Implementing material models
- [solvers_dev.md](solvers_dev.md) - Implementing solution algorithms
- [io_dev.md](io_dev.md) - Implementing I/O modules
- Available element types documentation
