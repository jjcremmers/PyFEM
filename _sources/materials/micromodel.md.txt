# MicroModel

The `MicroModel` material model implements a simple FE2-style constitutive response by solving a micro-scale representative volume element (RVE) inside each macro-scale material point. The macro model supplies the strain increment, the micro model solves one nonlinear RVE step with periodic boundary conditions, and the homogenized micro stress is returned as the macro stress.

## Overview
- **Material type:** `MicroModel`
- **Kinematics:** 2D small-strain only
- **Micro boundary conditions:** periodic rectangular RVE using the existing `RVE` model
- **Macro tangent:** numerical finite-difference tangent from repeated micro solves

## How It Works
At each integration point, `MicroModel`:

1. Reads a separate micro `.pro` input file in the constructor.
2. Creates one independent embedded micro problem per material point.
3. Applies the macro strain increment as `unitStrain` on the micro `RVE` model.
4. Runs one nonlinear micro solve.
5. Returns the homogenized stress from the micro model.
6. Approximates the consistent tangent numerically by perturbing the strain increment.

This makes the model easy to reuse with existing PyFEM elements, materials, solvers, and RVE constraints, but it is computationally expensive because every macro material evaluation triggers several micro analyses.

## Required Structure of the Micro Problem
The micro input file referenced by `rveFile` must:

- define `models = ["rve"]`
- contain an `rve` block with `type = "RVE"`
- define a rectangular RVE mesh
- provide the node groups `Left`, `Right`, `Top`, and `Bottom`

The current implementation is intended for 2D RVEs built from small-strain continuum elements.

## Parameters
### Mandatory
- `type`: Must be set to `"MicroModel"`
- `rveFile`: Path to the micro-scale `.pro` file

### Optional
- `perturbation`: Finite-difference perturbation used to compute the macro tangent. Default: `1.0e-7`
- `exportSamples`: List of micro-model sample IDs to export, a single integer, or `"all"`. Default: no export
- `exportFormat`: Export format for selected micro models. Accepts `"vtu"`, `"h5"`, or a list such as `["vtu","h5"]`. Default: `"vtu"`
- `exportPrefix`: Output basename for exported micro models. If omitted, a default name is generated from the micro input file

## Output
The material stores the homogenized in-plane stresses as nodal output fields:

- `S11`
- `S22`
- `S12`

## Default Export Naming
If `exportSamples` is set and `exportPrefix` is omitted, exported files are written next to the micro input file with the basename:

```text
<micro_rve_stem>_gp<sampleID>_macro_t<cycle>
```

For example, sample `0` at macro cycle `1` from `micro_rve.pro` becomes:

- `micro_rve_gp0_macro_t1.pvd`
- `micro_rve_gp0_macro_t1_t1.vtu`

If `exportFormat = "h5"`, the corresponding file is:

- `micro_rve_gp0_macro_t1.h5`

## Examples
The demonstrators in `examples/materials/micromodel/` include:

- `macro_fe2.pro`: minimal one-element macro test with a simple structured micro mesh and uniform extension
- `macro_bending.pro`: cantilever-style macro problem with a bending-dominated stress field and a two-fibre periodic composite RVE

### Minimal FE2 Example
The simplest setup uses one macro quadrilateral and one embedded micro RVE.

Macro input:
```text
MacroElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "MicroModel";
    rveFile = "examples/materials/micromodel/micro_rve.pro";
    perturbation = 1.0e-7;
    exportSamples = [0];
  };
};
```

Micro input:
```text
models = ["rve"];

rve =
{
  type = "RVE";
  boundaryType = "Periodic";
  unitStrain = [0.0,0.0,0.0];
};
```

Run the minimal example from the repository root:
```bash
python -m pyfem.core.cli examples/materials/micromodel/macro_fe2.pro
```

Related example files:
- `examples/materials/micromodel/macro_fe2.pro`
- `examples/materials/micromodel/macro_fe2.dat`
- `examples/materials/micromodel/micro_rve.pro`
- `examples/materials/micromodel/micro_rve.dat`

### Bending Example with Fibre Microstructure
The richer example uses a four-element macro cantilever and a periodic micro model built from a real fibre cross-section mesh with separate `Fibre` and `Epoxy` phases. Opposite horizontal tip forces create a bending moment, so the micro RVEs near the top and bottom of the beam are driven by clearly different macro strains.

Run it from the repository root:
```bash
python -m pyfem.core.cli examples/materials/micromodel/macro_bending.pro
```

Related example files:
- `examples/materials/micromodel/macro_bending.pro`
- `examples/materials/micromodel/macro_bending.dat`
- `examples/materials/micromodel/micro_rve_fibres.pro`
- `examples/materials/micromodel/micro_rve_fibres.dat`
- `examples/models/rve/composite/fibre.msh`

## Limitations
- Only 2D small-strain macro kinematics are supported.
- The micro problem must use the `RVE` model.
- The RVE is assumed to be rectangular with matching periodic boundary node groups.
- The macro tangent is numerical rather than analytical, so the model is relatively slow.
- Export selection is currently based on internal sample IDs, not on element numbers or physical coordinates.
