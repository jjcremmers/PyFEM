# PyFEM Installation Guide

PyFEM can be installed directly from the [GitHub source](https://github.com/jjcremmers/PyFEM).
Both the **Python API** and the **command-line interface (CLI)** are included.

## Requirements

**System Requirements:**
- Python 3.11 or newer
- [uv](https://docs.astral.sh/uv/) (recommended) or pip
- Git (for cloning the repository)

**Python Dependencies** (installed automatically):
- numpy
- scipy
- matplotlib
- meshio
- h5py
- PySide6
- vtk

## Installation with uv (Recommended)

[uv](https://docs.astral.sh/uv/) installs a compatible Python version, creates a virtual
environment, and installs PyFEM and its dependencies.

```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
uv sync
```

This creates a `.venv` directory and installs the `pyfem` and `pyfem-gui` commands.
Run them via `uv run`:

```bash
uv run pyfem --help
uv run pyfem-gui
```

To activate the environment manually:

```bash
source .venv/bin/activate  # Linux / macOS
.venv\Scripts\activate     # Windows
pyfem --help
```

### Development setup

`uv sync` also installs dev tools (pytest, coverage, ruff):

```bash
uv sync
uv run pytest
uv run coverage run -m pytest -q
uv run coverage report
uv run ruff check pyfem test
uv run ruff format --check pyfem test
uv build
```

## Installation with pip

If you prefer pip, create a virtual environment first:

```bash
python3 -m venv .venv
source .venv/bin/activate  # Linux / macOS
# .venv\Scripts\activate   # Windows
pip install .
```

For editable development installs:

```bash
pip install -e .
pip install pytest ruff
```

### Direct from GitHub

```bash
pip install git+https://github.com/jjcremmers/PyFEM.git
```

## Verifying the Installation

After installation, verify that both the CLI and API work correctly.

### Checking the CLI
```bash
pyfem --help
cd examples/ch02
pyfem PatchTest.pro
```
Expected output includes solver iterations, convergence information, and generated output files.

### Checking the API
```python
from pyfem import run
results = run("examples/ch02/PatchTest.pro")
print(results['globdat'].state)  # Displacement vector
```
If both tests complete without errors, the installation is successful.

## Using PyFEM

### Command-Line Interface (CLI)
**Basic Usage:**
```bash
pyfem input_file.pro
```
**Command-Line Options:**
```bash
pyfem --help                    # Show help
pyfem -i input.pro              # Specify input file
pyfem -d state.dump             # Restart from dump file
pyfem -p param=value            # Override parameter
```
**Examples:**
```bash
pyfem examples/ch03/cantilever8.pro
pyfem -d results_cycle100.dump
pyfem -i model.pro -p E=210000
```

### Python API
**Simple Usage - Run to Completion:**
```python
from pyfem import run
results = run("input.pro")
globdat = results['globdat']
displacements = globdat.state
props = results['props']
```
**Advanced Usage - Step-by-Step Control:**
```python
from pyfem.core.api import PyFEMAPI
api = PyFEMAPI("input.pro")
while api.isActive:
    api.step()
    current_disp = api.globdat.state
    load_factor = api.globdat.lam
    if load_factor > 5.0:
        print(f"Load factor reached {load_factor}")
results = api.getResults()
api.close()
```
**Loading from Pre-parsed Input:**
```python
from pyfem.io.InputReader import InputRead
from pyfem.core.api import PyFEMAPI
props, globdat = InputRead("input.pro")
props.solver.tol = 1e-6
api = PyFEMAPI((props, globdat))
api.runAll()
```
**Accessing Results:**
```python
globdat = api.globdat
displacements = globdat.state
node_coords = globdat.nodes.getNodeCoords(nodeID)
for name in globdat.outputNames:
    data = globdat.getData(name, node_list)
load_factor = globdat.lam
cycle = globdat.solverStatus.cycle
converged = globdat.solverStatus.converged
```

## Updating PyFEM
```bash
cd PyFEM
git pull origin main
uv sync
# Or with pip:
pip install --upgrade .
# Or if installed directly from GitHub
pip install --upgrade git+https://github.com/jjcremmers/PyFEM.git
```

## Uninstalling
```bash
pip uninstall pyfem
# Or remove the uv-managed environment:
rm -rf .venv
```

## Troubleshooting
**1. "pyfem: command not found"**
```bash
which pyfem  # Linux/macOS
where pyfem  # Windows
~/.local/bin/pyfem input.pro
# With uv, use:
uv run pyfem input.pro
```
**2. Import errors**
```bash
uv sync --reinstall
# Or with pip:
pip install --force-reinstall pyfem
```
**3. VTK, GUI, or documentation build issues**

```bash
sudo apt-get install -y libgl1 libxkbcommon-x11-0  # Linux (tests/GUI; libgl1 also needed for doc builds)
# On macOS, install XQuartz
brew install --cask xquartz
```
**4. Permission errors during installation**
Use a virtual environment (`uv sync` or `python -m venv .venv`) rather than installing system-wide.

## Platform-Specific Notes
**Linux:**
```bash
# uv installs its own Python; no system packages required for the venv
curl -LsSf https://astral.sh/uv/install.sh | sh
```
**macOS:**
```bash
brew install uv
brew install --cask xquartz  # for GUI / VTK
```
**Windows:**
1. Install [uv](https://docs.astral.sh/uv/getting-started/installation/)
2. Use PowerShell or Command Prompt
3. Install Git for Windows: [git-scm.com](https://git-scm.com/)

## Running Examples
```bash
cd examples
cd ch02
pyfem PatchTest.pro
cd ch03
pyfem cantilever8.pro
paraview cantilever8.pvd
```
Each example directory contains:
- `.pro` files: Input files
- `.dat` files: Mesh files
- Output files: VTK, text, plots

## Getting Help
- **Documentation**: https://pyfem.readthedocs.io/
- **GitHub Issues**: https://github.com/jjcremmers/PyFEM/issues
- **Examples**: See the `examples/` directory
- **Book**: "Non-Linear Finite Element Analysis of Solids and Structures" by de Borst et al., John Wiley & Sons, 2012

## Building documentation

Same install path as CI and Read the Docs:

```bash
uv sync --extra docs --no-dev
uv run sphinx-build -M html doc doc/_build
```

Open `doc/_build/html/index.html`. API reference is generated by `sphinx-autoapi` at build time (under `doc/api/`, gitignored).

On Linux, install `libgl1` if the build fails when loading VTK modules for viewcode.

## Next Steps
1. Read the [Quickstart guide](../introduction/quickstart.md)
2. Explore examples in the `examples/` directory
3. Review the [API reference](https://pyfem.readthedocs.io/en/latest/api/pyfem/index.html)
4. For development, see the [developer overview](../develop/overview.md)
