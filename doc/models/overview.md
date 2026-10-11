# Models

PyFEM includes specialized computational models that extend beyond standard finite element formulations. These models implement additional physics, constraints, or multi-scale approaches for specific analysis types.

## Introduction

Models in PyFEM provide:
- Periodic boundary conditions for representative volume elements (RVE)
- Contact mechanics for interface interactions
- Multi-scale coupling between different length scales
- Special constraint systems for specific applications

Models are configured in the `.pro` input file and run alongside the solver to apply constraints, manage coupling, or enforce special conditions.

Models are configured in the `.pro` input file using a model block:

```text
model_name = {
  type = "ModelType";
  # model-specific parameters
};
```

Multiple models can be active simultaneously, each handling different aspects
of the analysis.

## Contact

The Contact model implements penalty-based contact mechanics for interactions
between deformable bodies and rigid obstacles.

```{toctree}
:maxdepth: -1

contact.md
```

## RVE

The RVE model supports representative volume element analyses and periodic
boundary conditions for computational homogenization.

```{toctree}
:maxdepth: -1

rve.md
```

