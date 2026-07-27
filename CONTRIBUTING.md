# Contributing to PyFEM

Thank you for your interest in contributing to PyFEM.

PyFEM is an educational and research-oriented finite element code. Its main
purpose is to make finite element formulations, numerical algorithms, and
software structures clear and accessible. Contributions should therefore
prioritize correctness, readability, and educational value.

## Ways to Contribute

Contributions may include:

- reporting bugs;
- correcting or improving documentation;
- adding examples or exercises;
- improving tests;
- implementing new elements, material models, solvers, or output modules;
- improving code quality, usability, or portability.

Before starting a substantial change, please open a GitHub issue to describe
your proposal. This allows the maintainers and contributors to discuss the
scope and approach before significant work is invested.

## Reporting Bugs

Please use GitHub Issues to report bugs.

A useful bug report should include:

- a clear description of the problem;
- the PyFEM version or commit used;
- the operating system and Python version;
- the input files or a minimal example needed to reproduce the problem;
- the expected behaviour;
- the actual behaviour;
- relevant error messages or output.

Please remove confidential or proprietary data before attaching files.

Security vulnerabilities should not be reported in a public issue. Follow the
instructions in `SECURITY.md` instead.

## Requesting Features

Feature requests are welcome, particularly when they support education,
research, or the transparent implementation of finite element methods.

Please describe:

- the problem the feature would solve;
- its expected educational or scientific value;
- the proposed behaviour;
- relevant references, equations, or publications;
- possible alternatives, when applicable.

## Development Setup

Fork the repository on GitHub and clone your fork:

```bash
git clone https://github.com/<your-username>/PyFEM.git
cd PyFEM
```

Create and activate a virtual environment:

```bash
python -m venv pyfem-env
```

On Linux or macOS:

```bash
source pyfem-env/bin/activate
```

On Windows PowerShell:

```powershell
pyfem-env\Scripts\Activate.ps1
```

Install PyFEM in editable mode:

```bash
python -m pip install --upgrade pip
pip install -e .
```

Create a branch for your contribution:

```bash
git checkout -b feature/short-description
```

Use a focused branch name such as:

- `feature/add-new-element`;
- `fix/newton-convergence`;
- `docs/material-model-example`;
- `test/beam-elements`.

## Code Guidelines

PyFEM is intended to be readable by students, researchers, and developers.
Please keep implementations as clear and direct as reasonably possible.

Contributions should:

- follow PEP 8 where practical;
- use descriptive names for classes, functions, and variables;
- keep functions and classes focused on a clear responsibility;
- include docstrings for public classes and functions;
- explain non-obvious algorithms and finite element formulations;
- avoid unnecessary abstractions or external dependencies;
- remain consistent with the structure and conventions of the surrounding code.

Mathematical notation in the code and documentation should be consistent with
the notation used elsewhere in PyFEM. When implementing a published method,
please include an appropriate reference in the documentation or source code.

## Numerical Implementations

New finite element formulations and numerical algorithms should be accompanied
by enough information to assess their correctness.

Where applicable, include:

- the governing equations or formulation;
- assumptions and limitations;
- references to textbooks or publications;
- the adopted tensor or Voigt convention;
- units and sign conventions;
- details of the time-integration or solution procedure;
- a verification example with a known or independently established result.

Avoid introducing unexplained numerical constants. Parameters such as
tolerances, iteration limits, and stabilization factors should be documented.

## Tests and Verification

Contributions that change numerical behaviour should include suitable tests or
verification examples.

Depending on the contribution, this may include:

- unit tests for individual functions;
- regression tests for existing analyses;
- patch tests;
- comparison with analytical solutions;
- mesh-convergence studies;
- comparison with published benchmark results;
- checks of symmetry, conservation, or objectivity.

Run the available tests before submitting a pull request. If a test cannot be
added, explain in the pull request how the change was verified.

## Examples and Documentation

New user-facing functionality should include documentation and, where
appropriate, a small example.

Examples should:

- focus on one main concept;
- use reasonably small models;
- run in a practical amount of time;
- contain comments or accompanying explanations;
- avoid unnecessary output files;
- be suitable for educational use whenever possible.

Please update relevant documentation when changing an input format, command,
public API, element, material model, solver, or output module.

## Commit Guidelines

Keep commits focused and use clear commit messages.

Examples:

```text
Add plane-stress material test
Fix convergence check in Newton solver
Document periodic boundary conditions
Refactor VTK output without changing behaviour
```

Avoid combining unrelated changes in one commit. Do not commit generated
results, temporary files, editor settings, or large binary files unless they
are essential to the contribution.

## Pull Requests

Before opening a pull request:

1. Update your branch with the latest changes from `main`.
2. Run the relevant tests and examples.
3. Review your changes for clarity and unnecessary files.
4. Update documentation where needed.
5. Describe how the contribution was verified.

A pull request should explain:

- what was changed;
- why the change is useful;
- whether numerical results are affected;
- how the change was tested or verified;
- any limitations or remaining work;
- related GitHub issues.

Keep pull requests reasonably small and focused. Large changes may be easier to
review when divided into several coordinated pull requests.

## Review Process

Contributions are reviewed for:

- correctness;
- readability;
- educational value;
- consistency with the PyFEM architecture;
- test and documentation quality;
- impact on existing examples and users.

Review comments are intended to improve the contribution and maintain a clear,
reliable educational codebase. A contribution may require revisions before it
is accepted.

## Backward Compatibility

Please avoid breaking existing input files, examples, or public interfaces
without a strong reason.

When a breaking change is necessary:

- explain why it is needed;
- identify the affected functionality;
- update examples and documentation;
- provide migration guidance where practical.

## Licensing

By submitting a contribution, you agree that your contribution may be
distributed under the same license as PyFEM.

Only submit work that you have the right to contribute. Do not include
proprietary code, restricted data, or material that is incompatible with the
project license.

## Academic and AI-Assisted Contributions

Contributors remain responsible for the correctness, originality, licensing,
and documentation of all submitted work, including work created with the
assistance of generative AI tools.

Do not submit generated code that you do not understand or cannot verify.
AI-assisted contributions should be reviewed carefully for:

- numerical and logical correctness;
- fabricated functions, references, or APIs;
- licensing or attribution concerns;
- unnecessary complexity;
- inconsistencies with PyFEM conventions.

## Questions

For questions about using PyFEM, open a GitHub Discussion or, when appropriate,
a GitHub Issue.

For security-related matters, follow `SECURITY.md`.
