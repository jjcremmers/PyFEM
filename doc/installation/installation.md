# PyFEM Installation

PyFEM requires Python 3.9 or newer, `pip`, and Git. This page describes how to
install PyFEM from its source repository and provides platform-specific
instructions for Windows, Linux, WSL, and macOS.

The standard installation includes the core dependencies and VTK output
support. It also creates the `pyfem` command-line tool. A virtual environment
is recommended because it keeps PyFEM's dependencies separate from other
Python projects.

## Quick installation

### Windows (Command Prompt)

After installing Python and Git, open **Command Prompt** (`cmd.exe`) and run:

```cmd
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
py -m venv pyfem-env
pyfem-env\Scripts\activate.bat
python -m pip install --upgrade pip
python -m pip install .
```

If Python or Git is not already available, install Python 3.9 or newer from
[python.org](https://www.python.org/downloads/) and Git from
[git-scm.com](https://git-scm.com/). When installing Python, enable **Add
Python to PATH**.

### Linux and WSL

In a terminal, or in a WSL shell, run:

```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
python3 -m venv pyfem-env
source pyfem-env/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

On Debian or Ubuntu, install the system prerequisites first if necessary:

```bash
sudo apt-get update
sudo apt-get install python3 python3-venv python3-pip git
```

### macOS

In Terminal, run:

```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
python3 -m venv pyfem-env
source pyfem-env/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

If Python or Git is not already installed, install them with
[Homebrew](https://brew.sh/):

```bash
brew install python git
```

## Advanced installation

The virtual environment should be activated before running the commands in
this section. Use the `gui` extra to add graphical user interface support and
the `all` extra to install all optional dependencies.

### Windows (Command Prompt)

For an editable development installation, changes to the source are available
without reinstalling:

```cmd
cd PyFEM
pyfem-env\Scripts\activate.bat
python -m pip install -e ".[dev]"
```

For GUI support use `python -m pip install -e ".[dev,all]"`. To install the
latest source directly from GitHub instead:

```cmd
python -m pip install git+https://github.com/jjcremmers/PyFEM.git
```

### Linux and WSL

For an editable development installation, changes to the source are available
without reinstalling:

```bash
cd PyFEM
source pyfem-env/bin/activate
python -m pip install -e ".[dev]"
```

For GUI support use `python -m pip install -e ".[dev,all]"`. To install the
latest source directly from GitHub instead:

```bash
python -m pip install git+https://github.com/jjcremmers/PyFEM.git
```

On headless Linux systems, GUI or VTK issues may require these system libraries:

```bash
sudo apt-get install libgl1-mesa-glx libxkbcommon-x11-0
```

### macOS

For an editable development installation, changes to the source are available
without reinstalling:

```bash
cd PyFEM
source pyfem-env/bin/activate
python -m pip install -e ".[dev]"
```

For GUI support use `python -m pip install -e ".[dev,all]"`. To install the
latest source directly from GitHub instead:

```bash
python -m pip install git+https://github.com/jjcremmers/PyFEM.git
```

If graphical applications require it, install XQuartz:

```bash
brew install --cask xquartz
```

## Verify your installation

With the virtual environment active, run:

```bash
pyfem --version
pyfem --help
```

You can also run an example from the repository. From the repository root:

```bash
cd examples/ch02
pyfem PatchTest8.pro
```

To verify the Python API, run the following from the repository root:

```python
from pyfem import run

results = run("examples/ch02/PatchTest8.pro")
print(results["globdat"].state)
```

## Updating, uninstalling, and troubleshooting

From the cloned repository, update an existing installation with:

```bash
git pull origin main
python -m pip install --upgrade .
```

Uninstall PyFEM with:

```bash
python -m pip uninstall pyfem
```

If the `pyfem` command cannot be found, activate `pyfem-env` again. On Linux or
macOS, use `which pyfem` to show the executable location; in Windows Command
Prompt, use `where pyfem`. To repair the installation, run:

```bash
python -m pip install --force-reinstall .
```

For development, run the test suite with:

```bash
python -m pytest test/
```
