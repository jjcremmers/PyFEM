# OutputWriter

The `OutputWriter` I/O module writes a human-readable summary of all nodes per cycle to a text file. It uses the global print routine to include node IDs, coordinates, and active dof values.

## Overview
- **Module type:** `OutputWriter`
- Output file: `<prefix>_glob.out` by default, or a custom `filename`.
- On-screen print: optionally prints the same summary to stdout.

## Parameters
### Mandatory
- `type`: Must be set to `"OutputWriter"`

### Optional
(See documentation for details)
