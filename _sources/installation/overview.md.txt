# PyFEM Installation Guide

PyFEM can be installed directly from the [GitHub source](https://github.com/jjcremmers/PyFEM).
Both the **Python API** and the **command-line interface (CLI)** are included.

## Requirements

**System Requirements:**
- Python 3.9 or newer
- pip (Python package manager)
- Git (for cloning the repository)

**Python Dependencies** (installed automatically):
- numpy
- scipy
- matplotlib
- meshio
- h5py
- PySide6
- vtk

**Recommended: Virtual Environment**
It's recommended to install PyFEM in a virtual environment to avoid conflicts with other Python packages:

```bash
# Create virtual environment
python3 -m venv pyfem-env
# Activate on Linux / macOS
source pyfem-env/bin/activate
# Activate on Windows PowerShell
pyfem-env\Scripts\activate
# Activate on Windows Command Prompt
pyfem-env\Scripts\activate.bat
```

## Installation Steps

### Method 1: Standard Installation (Recommended)
```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
pip install .
```
This installs PyFEM and all dependencies, and creates the `pyfem` and `pyfem-gui` command-line executables.

### Method 2: Development Installation
```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
pip install -e .
```
The `-e` flag installs in "editable" mode, so changes to the source code are immediately reflected without reinstalling.

### Method 3: Direct from GitHub (Advanced)
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
~/.local/bin/pyfem input.pro
```
**2. Import errors**
```bash
pip install --force-reinstall pyfem
```
**3. VTK or GUI issues**
```bash
sudo apt-get install libgl1-mesa-glx libxkbcommon-x11-0  # Linux
# On macOS, install XQuartz
brew install --cask xquartz
```
**4. Permission errors during installation**
```bash
pip install --user .
```

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
1. Install Python 3.9+ from [python.org](https://www.python.org/downloads/)
2. Ensure "Add Python to PATH" is checked
3. Use PowerShell or Command Prompt
4. Install Git for Windows: [git-scm.com](https://git-scm.com/)

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

## Development Setup
```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
pip install -e .
python -m pytest test/
python -m black pyfem/
python -m mypy pyfem/
```

## Getting Help
- **Documentation**: https://pyfem.readthedocs.io/
- **GitHub Issues**: https://github.com/jjcremmers/PyFEM/issues
- **Examples**: See the `examples/` directory
- **Book**: "Non-Linear Finite Element Analysis of Solids and Structures" by de Borst et al., John Wiley & Sons, 2012

## Next Steps
1. Read the [Quickstart tutorial](../tutorials/quickstart.md)
2. Explore examples in the `examples/` directory
3. Review the [module documentation](../pyfem.md)
4. For development, see the [developer overview](../develop/overview.md)
