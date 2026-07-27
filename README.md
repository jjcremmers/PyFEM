# PyFEM: A Python Finite Element Code

[![Python Version](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![CI](https://github.com/jjcremmers/PyFEM/actions/workflows/ci.yml/badge.svg)](https://github.com/jjcremmers/PyFEM/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-brightgreen.svg)](https://pyfem.readthedocs.io/)
[![GitHub Stars](https://img.shields.io/github/stars/jjcremmers/PyFEM?style=social)](https://github.com/jjcremmers/PyFEM/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/jjcremmers/PyFEM)](https://github.com/jjcremmers/PyFEM/issues)
[![Cite](https://img.shields.io/badge/Cite-How%20to%20cite-blue.svg)](doc/introduction/introduction.md#how-to-cite)

PyFEM is a Python-based finite element code designed for educational and research purposes in computational solid mechanics. The code emphasizes clarity and readability, making it ideal for learning, teaching, and prototyping finite element methods for nonlinear analysis.

## ✨ Features

- **Comprehensive Element Library**: Continuum elements (2D/3D, small and finite strain), beam elements, plate elements, interface elements
- **Material Models**: Linear elastic, plasticity with hardening, cohesive zone models for fracture
- **Nonlinear Solvers**: Newton-Raphson, arc-length (Riks), explicit dynamics
- **I/O Capabilities**: VTK output for ParaView, HDF5, text-based formats
- **Python API**: Programmatic control for custom analyses and workflows
- **Multi-scale Modeling**: RVE (Representative Volume Element) with periodic boundary conditions
- **Well-Documented**: Extensive documentation with examples and developer guides

## 📚 About the Book

PyFEM accompanies the textbook:

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel<br>
[**Non-Linear Finite Element Analysis of Solids and Structures**](https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449)<br>
John Wiley and Sons, 2012, ISBN 978-0470666449

<img src="https://media.wiley.com/product_data/coverImage300/47/04706664/0470666447.jpg" width="200" alt="Book Cover">

The code is open source and intended for educational and scientific purposes. If you use PyFEM in your research, please cite the book in your publications.

## 📥 Download and Installation

### Requirements

- Python 3.11 or newer
- [uv](https://docs.astral.sh/uv/) (recommended) or pip
- Git (for cloning the repository)

### Quick Installation (recommended)

[uv](https://docs.astral.sh/uv/) manages Python, the virtual environment, and dependencies in one step:

```bash
# Clone the repository
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM

# Install Python 3.13 (if needed), create .venv, and install PyFEM
uv sync
```

Run PyFEM without activating the environment:

```bash
uv run pyfem --help
cd examples/ch02
uv run pyfem PatchTest.pro
```

### Development Installation

For contributors, `uv sync` installs PyFEM in editable mode together with dev tools (pytest, coverage, ruff):

```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
uv sync
uv run pytest
uv run coverage run -m pytest -q
uv run coverage report
uv run ruff check pyfem test
```

Build and verify the wheel (same as CI):

```bash
uv build
```

### Building documentation

CI and Read the Docs use `--no-dev` so only the `docs` extra is installed. Locally, omit `--no-dev` if you want dev tools in the same environment.

```bash
uv sync --extra docs --no-dev
uv run sphinx-build -M html doc doc/_build
```

On Linux, if the doc build fails loading VTK/OpenGL, install `libgl1` (see the [Installation Guide](doc/installation/overview.md)).

### Alternative: pip

If you prefer pip, use a virtual environment and install from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate  # Linux/macOS
# or: .venv\Scripts\activate  # Windows
pip install .
```

For detailed installation instructions including platform-specific notes, see the [Installation Guide](doc/installation/overview.md).

## 🚀 Quick Start

### Command-Line Interface

Run a PyFEM analysis from the command line:

```bash
# Navigate to examples
cd examples/ch02

# Run an example
uv run pyfem PatchTest.pro
```

View results in [ParaView](https://www.paraview.org/download/):

```bash
paraview PatchTest.pvd
```

## 📖 Documentation

### User Guide

- **[Installation Guide](doc/installation/overview.md)** - Complete installation instructions
- **[Quick Start Tutorial](doc/introduction/quickstart.md)** - Get started with PyFEM
- **[Elements](doc/elements/overview.md)** - Available element formulations
- **[Materials](doc/materials/overview.md)** - Material model documentation
- **[Solvers](doc/solvers/overview.md)** - Solution algorithms
- **[I/O Modules](doc/io/overview.md)** - Input/output capabilities
- **[Models](doc/models/overview.md)** - Special models (RVE, contact)
- **[Examples](examples/)** - Collection of example analyses

### Developer Guide

For contributors and those extending PyFEM:

- **[Developer's Overview](doc/develop/overview.md)** - Getting started with development
- **[Implementing Elements](doc/develop/elements_dev.md)** - Creating new element formulations
- **[Implementing Materials](doc/develop/materials_dev.md)** - Developing material models
- **[Implementing Solvers](doc/develop/solvers_dev.md)** - Creating solution algorithms
- **[Implementing I/O Modules](doc/develop/io_dev.md)** - Adding input/output capabilities

### API Reference

- **[API Documentation](https://pyfem.readthedocs.io/en/latest/api/pyfem/index.html)** - Python API reference (generated from source)

## 🎯 Example Gallery

PyFEM includes numerous examples organized by chapter from the book:

```bash
examples/
├── ch02/    # Linear elasticity and patch tests
├── ch03/    # Nonlinear analysis
├── ch04/    # Isoparametric elements
├── ch05/    # Element technology
├── ch06/    # Plasticity
├── ch09/    # Dynamics
├── ch13/    # Contact mechanics
├── ch15/    # Damage and fracture
├── elements/ # Element-specific examples
├── materials/ # Material model examples
├── models/   # Special models (RVE)
└── plate/    # Plate and shell examples
```

Each directory contains input files (`.pro`), mesh files (`.dat`), and generates output files for visualization.

## 🤝 Contributing

Contributions are welcome! Please see the [Developer's Guide](doc/develop/overview.md) for:

- Code style and conventions
- Testing guidelines
- Documentation requirements
- Pull request process

To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes and add tests
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

## 📝 License

PyFEM is distributed under the MIT License. See [LICENSE](LICENSE) for details.

## 📧 Contact and Support

- **Issues**: [GitHub Issues](https://github.com/jjcremmers/PyFEM/issues)
- **Discussions**: [GitHub Discussions](https://github.com/jjcremmers/PyFEM/discussions)
- **Email**: Contact through GitHub

## 🌟 Citing PyFEM

If PyFEM contributes to a publication, please cite:

**Software:**
```
J.J.C. Remmers (2026). PyFEM - A Python Finite Element Code.
https://github.com/jjcremmers/PyFEM
```

**Textbook:**
```
R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel (2012).
Non-Linear Finite Element Analysis of Solids and Structures, 2nd Edition.
John Wiley & Sons, ISBN 978-0470666449.
```

## 🙏 Acknowledgments

PyFEM is developed and maintained by:

- **Joris Remmers** - Eindhoven University of Technology
- Contributors and users from the computational mechanics community

The code accompanies the textbook by de Borst, Crisfield, Remmers, and Verhoosel, which provides the theoretical foundation for the implemented methods.

[paraViewURL]: paraview.org