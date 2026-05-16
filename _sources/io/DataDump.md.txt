# DataDump

The `DataDump` I/O module serializes the entire analysis state (`props` and `globdat`) to a pickle file for later restart or inspection.

## Overview
- **Module type:** `DataDump`
- Output file: `<prefix>_<cycle>.dump` or `<prefix>.dump` when `lastOnly`.
- Contents: both the properties object and global data.

## Parameters
### Mandatory
- `type`: Must be set to `"DataDump"`

### Optional
(See documentation for details)
