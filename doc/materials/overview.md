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

These models describe elasticity, plasticity, anisotropy, and multiscale
constitutive response in continuum materials.

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
PlaneStrainDamage.md
sandwichcore.md
transverseisotropic.md
ViscoElasticity.md
ViscoPlasticity.md
SOVS.md
```

## Cohesive Zone Models

Cohesive constitutive models describe the local relation between the
displacement jump across two adjoining surfaces and the tractions transferred
between them. The response generally includes an initial elastic or penalty
regime, followed by damage and softening as the interface opens or slides, and
may end in complete decohesion. Parameters such as the peak traction and
fracture energy control the onset and evolution of failure. These models are
used for delamination, adhesive failure, peeling, and other debonding
processes. They must be assigned within the material block of the
[`Interface` element](../elements/interface.md), which supplies the interface
kinematics and transfers the resulting normal and tangential tractions to the
adjacent continuum parts.

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
MaxStrain.md
MaxStress.md
TsaiWu.md
vonmises.md
```

