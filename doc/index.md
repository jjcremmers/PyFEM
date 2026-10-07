# PyFEM: A Python Finite Element Code

PyFEM is a Python-based finite element code designed for educational and research purposes in computational solid mechanics. The code emphasizes clarity and readability, making it ideal for learning, teaching, and prototyping finite element methods for nonlinear analysis.

## Features

- **Comprehensive Element Library**: Continuum elements (2D/3D, small and finite strain), beam elements, plate elements, interface elements
- **Material Models**: Linear elastic, plasticity with hardening, cohesive zone models for fracture
- **Nonlinear Solvers**: Newton-Raphson, arc-length (Riks), explicit dynamics
- **I/O Capabilities**: VTK output for ParaView, HDF5, text-based formats
- **Python API**: Programmatic control for custom analyses and workflows
- **Multi-scale Modeling**: RVE (Representative Volume Element) with periodic boundary conditions
- **Well-Documented**: Extensive documentation with examples and developer guides

## About the Book

PyFEM accompanies the textbook:

R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel<br>
[**Non-Linear Finite Element Analysis of Solids and Structures**](https://www.wiley.com/en-us/Nonlinear+Finite+Element+Analysis+of+Solids+and+Structures%2C+2nd+Edition-p-9780470666449)<br>
John Wiley and Sons, 2012, ISBN 978-0470666449

<img src="https://media.wiley.com/product_data/coverImage300/47/04706664/0470666447.jpg" width="200" alt="Book Cover">

The code is open source and intended for educational and scientific purposes. If you use PyFEM in your research, please cite the book in your publications.

## Download and Installation

The following is required for installation:

- Python 3.9 or higher
- pip package manager
- Git (for cloning the repository)

Quick Installation

```bash
# Clone the repository
git clone https://github.com/jjcremmers/PyFEM.git
cd PyFEM

# Create and activate a virtual environment on Linux/macOS
python3 -m venv pyfem-env
source pyfem-env/bin/activate

# On Windows (Command Prompt), use:
# py -m venv pyfem-env
# pyfem-env\\Scripts\\activate.bat

# Install PyFEM with pip
pip install .
```

The base installation includes VTK output support for ParaView. The graphical
user interface dependency is optional.

For detailed installation instructions including platform-specific notes
and development releases, see the [Installation Guide](installation/overview.md).

## Quick Start

Run a PyFEM analysis from the command line:

```bash
# Navigate to examples
cd examples/ch02

# Run an example
pyfem PatchTest8.pro
```

Useful command-line options:

```bash
pyfem input.pro                 # specify input file
pyfem --help                    # Show help
pyfem --version                 # Show installed PyFEM version
pyfem -d state.dump             # Restart from dump file
pyfem input.pro -p param=value  # Override parameter in input file
pyfem -test                     # Run the unit tests
pyfem -coverage                 # Run tests and print coverage
```

If you use the VTK writer to write the output to the harddisk, you can 
view the results using [ParaView](https://www.paraview.org/download/):

```bash
paraview input.pvd
```

## Example Gallery

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
└── models/   # Special models (RVE)
```

Each directory contains input files (`.pro`), mesh files (`.dat`), and generates output files for visualization.

## Contributing

Contributions are welcome. Please read
[CONTRIBUTING](CONTRIBUTING.md) before reporting a bug, proposing a
feature, or submitting a pull request.

## License

PyFEM is distributed under the MIT License. See [LICENSE](../LICENSE) for details.

## Contact and Support

- **Issues**: [GitHub Issues](https://github.com/jjcremmers/PyFEM/issues)
- **Discussions**: [GitHub Discussions](https://github.com/jjcremmers/PyFEM/discussions)
- **Email**: Contact through GitHub

## Citing PyFEM

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

## Acknowledgments

PyFEM is developed and maintained by:

- **Joris Remmers** - Eindhoven University of Technology
- Contributors and users from the computational mechanics community

The code accompanies the textbook by de Borst, Crisfield, Remmers, and Verhoosel, which provides the theoretical foundation for the implemented methods.

[paraViewURL]: paraview.org

```{toctree}
:maxdepth: 1

home.md
introduction/introduction.md
introduction/quickstart.md
installation/overview.md
usermanual.md
develop/overview.md
api.md
```
