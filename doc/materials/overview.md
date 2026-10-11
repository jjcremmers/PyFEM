# Materials

PyFEM provides a range of material models for linear and nonlinear analysis, including elastic, elastoplastic, and cohesive zone models. Material models define the constitutive behavior that relates stress to strain (or traction to displacement jump for cohesive elements).

## Introduction

Material models are configured within element definitions in the `.pro` input file. Each element group can have its own material definition, or materials can be defined globally and referenced by name.

Basic Syntax:
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

## Continuum Constitutive Models

Continuum constitutive models define how bulk materials respond to deformation.
They relate strain to stress and may represent elastic, plastic, anisotropic,
or multiscale behavior in solid continuum elements.

```{toctree}
:maxdepth: -1

crystal.md
hookeslaw.md
isokinematichardening.md
isotropic.md
isotropichardeningplasticity.md
micromodel.md
multimaterial.md
planestrain.md
planestress.md
sandwichcore.md
transverseisotropic.md
```

## Cohesive Zone Models

These models define traction–separation laws for interfaces, delamination,
debonding, and other cohesive-zone applications.

```{toctree}
:maxdepth: -1

dummy.md
powerlawmodei.md
thoulessmodei.md
xuneedleman.md
```

## Failure Models

These models evaluate stress-based or composite-specific failure criteria.

```{toctree}
:maxdepth: -1

larc03.md
vonmises.md
```

