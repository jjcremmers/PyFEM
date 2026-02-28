# Materials

PyFEM provides a range of material models for linear and nonlinear analysis, including elastic, elastoplastic, and cohesive zone models. Material models define the constitutive behavior that relates stress to strain (or traction to displacement jump for cohesive elements).

## Configuration in Input Files
Material models are configured within element definitions in the `.pro` input file. Each element group can have its own material definition, or materials can be defined globally and referenced by name.

### Basic Syntax
The material block is nested inside an element definition:
```text
ElementGroup = {
  type = "SmallStrainContinuum";
  material = {
    type = "PlaneStress";
    E    = 100.0;
    nu   = 0.3;
  };
};
```
The `type` parameter specifies which material model to use, and subsequent parameters define the material properties (e.g., Young's modulus `E`, Poisson's ratio `nu`).
