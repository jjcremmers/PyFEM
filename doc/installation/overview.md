# PyFEM Installation Guide

PyFEM can be installed directly from the [GitHub source](https://github.com/jjcremmers/PyFEM).
Both the **Python API** and the **command-line interface (CLI)** are included.

## Requirements

**System Requirements:**
- Python 3.10 or newer
- [uv](https://docs.astral.sh/uv/) for the recommended workflow, or pip
- Git (for cloning the repository)

**Python Dependencies** (installed automatically):
- numpy
- scipy
- matplotlib
- meshio
- h5py
- vtk

**Optional Dependencies:**
- `PySide6`: graphical user interface support

## Installation with uv

From the repository root, uv creates a virtual environment and installs PyFEM:

```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
uv sync
uv run pyfem --help
```

The package is installed in editable mode. To add the optional GUI:

```bash
uv sync --extra gui
uv run --no-sync pyfem-gui
```

### Development setup

```bash
uv sync --extra dev
uv run --no-sync pytest -q
uv build
```

Add `--extra gui` to `uv sync` when developing the GUI.

## Installation with pip

Create a virtual environment first. On Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install .
```

On Windows Command Prompt:

```text
py -m venv .venv
.venv\Scripts\activate.bat
python -m pip install .
```

This installs the core package and the `pyfem` command. When using pip in an
activated environment, run `pyfem` directly. For the GUI, run:

```bash
pip install ".[gui]"
```

For an editable development installation with test tools:

```bash
pip install -e ".[dev]"
```

### Direct from GitHub
```bash
pip install git+https://github.com/jjcremmers/PyFEM.git
```

## Verifying the Installation

After installation, verify that both the CLI and API work correctly.

### Checking the CLI
```bash
uv run pyfem --help
cd examples/ch02
uv run pyfem PatchTest8.pro
```
Expected output includes solver iterations, convergence information, and generated output files.

### Checking the API
```python
from pyfem import run
results = run("examples/ch02/PatchTest8.pro")
print(results['globdat'].state)  # Displacement vector
```
If both tests complete without errors, the installation is successful.

## Using PyFEM

### Command-Line Interface (CLI)
**Basic Usage:**
```bash
uv run pyfem input_file.pro
```
**Command-Line Options:**
```bash
uv run pyfem --help                    # Show help
uv run pyfem --version                 # Show installed PyFEM version
uv run pyfem -i input.pro              # Specify input file
uv run pyfem -d state.dump             # Restart from dump file
uv run pyfem -p param=value            # Override parameter
```
**Examples:**
```bash
cd examples/ch03
uv run pyfem cantilever8.pro
uv run pyfem -d results_cycle100.dump
uv run pyfem -i model.pro -p E=210000
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
while api.is_active:
    api.step()
    current_disp = api.globdat.state
    load_factor = api.globdat.lam
    if load_factor > 5.0:
        print(f"Load factor reached {load_factor}")
results = api.get_results()
api.close()
```
**Loading from Pre-parsed Input:**
```python
from pyfem.io.InputReader import InputRead
from pyfem.core.api import PyFEMAPI
props, globdat = InputRead("input.pro")
props.solver.tol = 1e-6
api = PyFEMAPI((props, globdat))
api.run_all()
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
# Or, with pip in an activated environment:
pip install --upgrade .
# Or if installed directly from GitHub
pip install --upgrade git+https://github.com/jjcremmers/PyFEM.git
```

## Uninstalling
```bash
pip uninstall pyfem
```

## Troubleshooting
**1. "pyfem: command not found"**
```bash
which pyfem  # Linux/macOS
where pyfem  # Windows
uv run pyfem input.pro
```
**2. Import errors**
```bash
uv sync --reinstall
# Or, with pip from the repository root:
pip install --force-reinstall .
```
**3. VTK or GUI issues**
```bash
sudo apt-get install libgl1 libxkbcommon-x11-0  # Linux
# On macOS, install XQuartz
brew install --cask xquartz
```
**4. Permission errors during installation**
Use `uv sync` or an activated virtual environment instead of installing into
the system Python.

## Platform-Specific Notes
**Linux:**
```bash
sudo apt-get install python3-venv python3-pip  # Debian/Ubuntu
sudo dnf install python3-virtualenv python3-pip  # Fedora/RHEL
```
**macOS:**
```bash
brew install python@3.11
brew install --cask xquartz
```
**Windows:**
1. Install Python 3.10+ from [python.org](https://www.python.org/downloads/)
2. Ensure "Add Python to PATH" is checked
3. Open **Command Prompt** (`cmd.exe`)
4. Install Git for Windows: [git-scm.com](https://git-scm.com/)
5. Run the installation commands from Command Prompt:

   ```text
   git clone https://github.com/jjcremmers/PyFEM.git
   cd PyFEM
   py -m venv pyfem-env
   pyfem-env\Scripts\activate.bat
   python -m pip install .
   ```

## Running Examples
```bash
cd examples
cd ch02
uv run pyfem PatchTest8.pro
cd ../ch03
uv run pyfem cantilever8.pro
paraview cantilever8.pvd
```
Each example directory contains:
- `.pro` files: Input files
- `.dat` files: Mesh files
- Output files: VTK, text, plots

## Building documentation

```bash
uv sync --extra docs
uv run --no-sync sphinx-build -b html doc doc/_build/html
```

Open `doc/_build/html/index.html`. The API reference is generated from source
by Sphinx. On Linux, install `libgl1` if VTK cannot load during the build.

## Getting Help
- **Documentation**: https://pyfem.readthedocs.io/
- **GitHub Issues**: https://github.com/jjcremmers/PyFEM/issues
- **Examples**: See the `examples/` directory
- **Book**: "Non-Linear Finite Element Analysis of Solids and Structures" by de Borst et al., John Wiley & Sons, 2012

## Next Steps
1. Read the [Quickstart guide](../introduction/quickstart.md)
2. Explore examples in the `examples/` directory
3. Review the [generated API reference](https://pyfem.readthedocs.io/en/latest/api/pyfem/index.html)
4. For development, see the [developer overview](../develop/overview.md)
