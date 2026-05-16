# InputReader

`InputReader` reads a `.pro` input file and constructs the analysis state (`props` and `globdat`). It can also read a previously saved dump to restore the full state.

## Overview
- **Function:** `InputReader(argv)` and `InputRead(fname, dname=None, parameters=None)`
- Reads properties from a `.pro` file; parses mesh, elements, dofs, models.
- Sets up logging, prefix, and initial global state.
- Optionally reads a dump file (pickle) to restore state instead of parsing.

## Command-Line Arguments
- `-i`, `--input`: Path to the `.pro` file
- `-d`, `--dump`: Path to a pickle dump (created by `DataDump`) to restore state
- `-p`, `--param`: Override parameters as `name=value` (can be repeated)
- `-h`, `--help`: Show help

## Programmatic Use
(See documentation for details)
