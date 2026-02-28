# GraphWriter

The `GraphWriter` I/O module writes column-based data per cycle and can optionally plot the first two columns live to a PNG. It supports data from `globdat.outputNames`, attributes on `globdat` (including arrays indexed by node/dof), and solver status fields.

## Overview
- **Module type:** `GraphWriter`
- Output file: `<prefix>.out` by default, or a custom `filename`.
- Live plot: when `onScreen = true`, saves a PNG `<prefix>.png` each update.
- Columns: configured via a list of column names; each column may have sub-parameters.

## Parameters
### Mandatory
- `type`: Must be set to `"GraphWriter"`

### Optional
(See documentation for details)
