# PyFEM: A Python Finite Element Code

[![Python Version](https://img.shields.io/badge/python-3.9--3.14-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://github.com/jjcremmers/PyFEM/tree/main/doc)
[![GitHub Stars](https://img.shields.io/github/stars/jjcremmers/PyFEM?style=social)](https://github.com/jjcremmers/PyFEM/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/jjcremmers/PyFEM)](https://github.com/jjcremmers/PyFEM/issues)
[![Cite](https://img.shields.io/badge/Cite-How%20to%20cite-blue.svg)](CITATION.cff)

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

For a quick setup, see the [Quick installation instructions](doc/installation/installation.md#quick-installation).
For all platform-specific instructions and advanced development setup, see the
[Installation Guide](doc/installation/installation.md).

## 🚀 Quick Start

For command-line usage, examples, input-file basics, visualization, and test
shortcuts, see the [Quickstart guide](doc/introduction/quickstart.md).

## 📖 Documentation

### User Guide

- **[Installation Guide](doc/installation/installation.md)** - Complete installation instructions
- **[Quickstart Guide](doc/introduction/quickstart.md)** - Run your first example
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

- **[API Documentation](doc/introduction/api.md)** - Python API reference
- **[Documentation Index](doc/index.md)** - Full documentation table of contents

## 🤝 Contributing

Contributions are welcome. Please read
[CONTRIBUTING.md](CONTRIBUTING.md) before reporting a bug, proposing a
feature, or submitting a pull request.

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
