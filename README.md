# PyFEM: A Python Finite Element Code

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://github.com/jjcremmers/PyFEM/tree/main/doc)
[![GitHub Stars](https://img.shields.io/github/stars/jjcremmers/PyFEM?style=social)](https://github.com/jjcremmers/PyFEM/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/jjcremmers/PyFEM)](https://github.com/jjcremmers/PyFEM/issues)
[![Cite](https://img.shields.io/badge/Cite-How%20to%20cite-blue.svg)](doc/index.rst#how-to-cite)

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

- Python 3.9 or higher
- pip package manager
- Git (for cloning the repository)

### Quick Installation

```bash
# Clone the repository
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM

# Install with pip
pip install .
```

### Development Installation

For developers who want to make changes and test immediately:

```bash
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM
pip install -e .
```

### Virtual Environment (Recommended)

```bash
# Create and activate virtual environment
python3 -m venv pyfem-env
source pyfem-env/bin/activate  # Linux/macOS
# or
pyfem-env\Scripts\activate  # Windows

# Install PyFEM
pip install .
```

For detailed installation instructions including platform-specific notes, see the [Installation Guide](doc/installation/overview.rst).

## 🚀 Quick Start

### Command-Line Interface

Run a PyFEM analysis from the command line:

```bash
# Navigate to examples
cd examples/ch02

# Run an example
pyfem PatchTest.pro
```

View results in [ParaView](https://www.paraview.org/download/):

```bash
paraview PatchTest.pvd
```

## 📖 Documentation

### User Guide

- **[Installation Guide](doc/installation/overview.rst)** - Complete installation instructions
- **[Quick Start Tutorial](doc/tutorials/quickstart.rst)** - Get started with PyFEM
- **[Elements](doc/elements/overview.rst)** - Available element formulations
- **[Materials](doc/materials/overview.rst)** - Material model documentation
- **[Solvers](doc/solvers/overview.rst)** - Solution algorithms
- **[I/O Modules](doc/io/overview.rst)** - Input/output capabilities
- **[Models](doc/models/overview.rst)** - Special models (RVE, contact)
- **[Examples](examples/)** - Collection of example analyses

### Developer Guide

For contributors and those extending PyFEM:

- **[Developer's Overview](doc/develop/overview.rst)** - Getting started with development
- **[Implementing Elements](doc/develop/elements_dev.rst)** - Creating new element formulations
- **[Implementing Materials](doc/develop/materials_dev.rst)** - Developing material models
- **[Implementing Solvers](doc/develop/solvers_dev.rst)** - Creating solution algorithms
- **[Implementing I/O Modules](doc/develop/io_dev.rst)** - Adding input/output capabilities

### API Reference

- **[API Documentation](doc/api.rst)** - Python API reference
- **[Module Documentation](doc/modules.rst)** - Complete module documentation

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

## 🛠️ Key Capabilities

### Nonlinear Analysis

- Geometrically nonlinear formulations
- Material nonlinearity (plasticity, damage)
- Path-following methods (arc-length)
- Contact mechanics

### Dynamic Analysis

- Explicit time integration
- Implicit dynamics
- Modal analysis
- Eigenvalue problems

### Multi-Scale Modeling

- Representative Volume Elements (RVE)
- Periodic boundary conditions
- Computational homogenization

### Fracture Mechanics

- Cohesive zone models
- Interface elements
- Crack propagation

## 🤝 Contributing

Contributions are welcome! Please see the [Developer's Guide](doc/develop/overview.rst) for:

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