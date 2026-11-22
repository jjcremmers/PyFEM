# Materials

PyFEM provides material models for linear and nonlinear analysis, including
elastic, elastoplastic, anisotropic, cohesive-zone, and failure models. A
material model defines the constitutive behavior that relates stress to strain,
or traction to displacement jump for cohesive elements. Materials are
configured within element definitions in the `.pro` input file. Each element
group can have its own material definition, or materials can be defined
globally and referenced by name.

## Basic syntax

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

The `type` parameter selects the material model. The remaining parameters
define its material properties, such as Young's modulus `E` and Poisson's
ratio `nu`. Refer to the individual material pages for required parameters,
units, assumptions, and examples.

## Continuum Constitutive Models

These models describe the constitutive response of bulk materials under
elastic, plastic, anisotropic, thermal, and multiscale assumptions. They cover
linear elasticity, isotropic and kinematic hardening, finite-strain behavior,
damage, viscoelasticity, viscoplasticity, composite and sandwich materials,
and specialized crystal or micromechanical models. Select a model according
to the material behavior, strain regime, internal variables, and element
formulation used in the analysis. The individual pages document the supported
parameters and their physical interpretation.

```{toctree}
:maxdepth: 1

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
:maxdepth: 1

dummy.md
powerlawmodei.md
thoulessmodei.md
xuneedleman.md
```

## Failure Models

These models evaluate stress-based or composite-specific failure criteria.
They identify the onset of failure by comparing computed stresses or strains
with material limits, and can distinguish common modes such as tension,
compression, shear, or interaction between them. The collection includes
general stress and strain limits, von Mises yielding, and composite criteria
such as Tsai–Wu and LaRC03. Failure models are used to assess a simulation's
state and should be selected consistently with the material orientation,
loading conditions, and intended failure definition.

```{toctree}
:maxdepth: 1

larc03.md
MaxStrain.md
MaxStress.md
TsaiWu.md
vonmises.md
```

