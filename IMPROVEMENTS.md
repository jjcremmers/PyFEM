# PyFEM Repository Improvement Backlog

This document collects suggested improvements from a repository audit focused
on installation, packaging, documentation, testing, continuous integration, and
repository hygiene. Items are grouped by priority and include concrete next
steps.

## High Priority

### 1. Fix the public Python API naming mismatch

The programmatic API currently mixes camelCase and snake_case names. The class
defines `isActive`, `runAll`, and `getResults`, while internal calls and
docstrings refer to `is_active`, `run_all`, and `get_results`. As a result,
`from pyfem import run` imports successfully but fails when called.

Recommended changes:

- Choose one naming style for the public API. For a Python package, prefer
  snake_case: `is_active`, `run_all`, and `get_results`.
- Keep backward-compatible aliases for existing users:
  `isActive`, `runAll`, and `getResults`.
- Add a small API regression test that calls `pyfem.run(...)` on a minimal
  example and checks that a result dictionary is returned.
- Update all documentation examples to use the chosen API names.

### 2. Repair README documentation links

The README links point to `.rst` files, but the documentation files in `doc/`
are Markdown files. Several links therefore lead to missing files on GitHub.

Recommended changes:

- Replace README links from `.rst` to `.md`, for example
  `doc/installation/overview.md`.
- Remove or replace links to missing files such as `doc/api.rst` and
  `doc/modules.rst`.
- Link the documentation badge to the published documentation site if one is
  available, not to the raw `doc/` directory.
- Add a link-checking job in CI, or run a documentation link checker locally
  before releases.

### 3. Make license metadata consistent

The repository license file and README state MIT, while package metadata in
`setup.py` declares GPL-3.0-only.

Recommended changes:

- Decide the intended project license.
- Make `LICENSE`, README badge/text, package metadata, classifiers, and SPDX
  source headers agree.
- If the project is MIT, update `setup.py` to `license="MIT"` and use the MIT
  classifier.
- If the project is GPL, replace the `LICENSE` file and README badge/text.
- Add a short note in `CONTRIBUTING.md` that contributions are accepted under
  the project license.

### 4. Use one source of truth for the package version

The package metadata declares version `3.00`, while `pyfem.__version__` is
`0.1.0`.

Recommended changes:

- Adopt a single canonical version source.
- Prefer modern package metadata in `pyproject.toml`.
- If keeping `setup.py`, read the version from one module or generate
  `pyfem.__version__` from installed metadata.
- Use a normalized version such as `3.0.0` instead of `3.00`.
- Add a test that checks `pyfem.__version__` matches installed package
  metadata.

### 5. Correct the quick-start example

The README asks users to run `examples/ch02/PatchTest.pro`, but that file is
not present. `examples/ch02/PatchTest8.pro` exists and runs successfully.

Recommended changes:

- Update the README quick-start command to use an existing example.
- Prefer a small, fast, robust example such as `examples/ch02/PatchTest8.pro`.
- Add a CI smoke test that runs the documented quick-start example.
- Keep generated output from that smoke test out of version control.

## Packaging and Installation

### 6. Add `pyproject.toml`

The project currently relies on legacy `setup.py` metadata.

Recommended changes:

- Add `pyproject.toml` with build-system metadata using `setuptools`.
- Move project metadata, dependencies, optional dependencies, package data, and
  console scripts into `pyproject.toml`.
- Keep a minimal compatibility `setup.py` only if needed.
- Configure common tools in `pyproject.toml`, such as pytest and coverage.

Suggested optional dependency groups:

- `pyfem[gui]`: `PySide6`
- `pyfem[vtk]`: `vtk`
- `pyfem[docs]`: Sphinx, MyST, theme, and documentation dependencies
- `pyfem[dev]`: test, coverage, formatting, linting, and type-checking tools

### 7. Make heavy GUI and visualization dependencies optional

Every installation currently requires `PySide6` and `vtk`. These are large
dependencies and can cause installation problems on servers, CI systems, and
headless environments.

Recommended changes:

- Move GUI-only dependencies to a `gui` extra.
- Move VTK-specific output support to a `vtk` extra if the core solver can run
  without it.
- Ensure imports of optional dependencies happen lazily, inside the modules or
  code paths that need them.
- Provide clear error messages when a user requests GUI or VTK functionality
  without the required extra installed.
- Document installation commands such as `pip install ".[gui,vtk]"`.

### 8. Fix Read the Docs dependency configuration

`.readthedocs.yaml` requests a `docs` extra, but the package does not define
one. The documentation requirements also omit `myst_parser`, even though
`doc/conf.py` uses it.

Recommended changes:

- Define a `docs` extra in package metadata, or remove `extra_requirements:
  docs` from `.readthedocs.yaml`.
- Add `myst_parser` and `sphinx_rtd_theme` to the canonical docs dependency
  list.
- Keep GitHub Actions documentation dependencies and Read the Docs dependencies
  aligned.
- Make documentation builds fail on warnings once the current warnings are
  cleaned up.

## Testing and Continuous Integration

### 9. Make pytest work

`python -m unittest discover -s test -p "*.py"` runs the suite, but
`python -m pytest test` collects zero tests. Many contributors will try pytest
first and get a false signal.

Recommended changes:

- Either rename test files/classes/methods to pytest-compatible discovery
  conventions or configure pytest explicitly.
- Add a `pyproject.toml` section such as:

```toml
[tool.pytest.ini_options]
python_files = ["test*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
testpaths = ["test"]
```

- Update `CONTRIBUTING.md` to document the canonical test command.
- Prefer one test runner in CI to avoid split behavior.

### 10. Test the installed package in CI

The current CI manually installs dependencies and runs tests from the source
tree. It does not verify that package metadata, dependencies, or console script
entry points work.

Recommended changes:

- Install the package in editable mode in CI with `python -m pip install -e
  ".[dev]"`.
- Add a packaging smoke test:
  `python -c "import pyfem; print(pyfem.__version__)"`.
- Add console-script smoke tests:
  `pyfem --help` and `pyfem examples/ch02/PatchTest8.pro`.
- Run CI on pull requests and pushes to active development branches.
- Consider adding a wheel build check with `python -m build`.

### 11. Align CI dependencies with package metadata

The CI dependency list is hand-written and differs from package metadata.

Recommended changes:

- Install dependencies through the package itself instead of duplicating the
  list in the workflow.
- If optional extras are added, test at least:
  core install, docs install, and full install.
- Add dependency caching only after the install path is stable.

## Command-Line Interface

### 12. Replace custom `getopt` parsing with `argparse`

The CLI advertises options such as `-i`, `-d`, and `-p`, but the help output
only shows a positional input file. The current parser also skips elements of
`argv` in a way that is easy to misuse.

Recommended changes:

- Use `argparse` in `pyfem.core.cli`.
- Support both `pyfem input.pro` and `pyfem -i input.pro`.
- Add documented options for restart dumps and parameter overrides.
- Make `pyfem --version` print the package version.
- Return clear user-facing errors for missing input files.
- Add CLI tests for positional input, `-i`, `--help`, and invalid input.

## Documentation Quality

### 13. Keep README, user manual, and CLI help synchronized

The README, installation guide, and CLI help currently describe overlapping but
not identical usage.

Recommended changes:

- Treat the README as a short entry point: install, one working example, docs
  link, citation, license.
- Put detailed CLI/API documentation in `doc/`.
- Generate or test CLI help examples where practical.
- Add a "tested with" section for supported Python versions and operating
  systems.
- Add troubleshooting notes for common optional dependency issues, especially
  VTK and GUI installation.

### 14. Add documentation build checks locally and in CI

The documentation setup uses Sphinx with Markdown through MyST.

Recommended changes:

- Add a documented command:
  `python -m sphinx -W -b html doc doc/_build/html`.
- Add missing docs dependencies to one canonical location.
- Gradually enable warnings-as-errors.
- Consider using `sphinx.ext.intersphinx` for Python, NumPy, SciPy, and
  Matplotlib references.
- Keep generated API documentation out of source directories unless
  intentionally committed.

## Repository Hygiene

### 15. Remove generated analysis outputs from version control

Many generated `.vtu`, `.pvd`, `.h5`, and `*_glob.out` files are present under
`examples/` and `exercises/`, even though `.gitignore` already ignores several
of these patterns.

Recommended changes:

- Decide whether generated outputs are reference artifacts or disposable
  outputs.
- If disposable, remove them from version control with `git rm --cached`.
- If some are needed as reference outputs, move them to a clearly named
  regression-data directory and document their purpose.
- Extend `.gitignore` to cover all generated output patterns, including
  generated `.h5` files if appropriate.
- Add tests that regenerate selected outputs and compare important scalar
  values rather than committing every visualization file.

### 16. Remove accidental backup files

Files such as `pyfem/io/MeshWriter.py1` and `pyfem/util/vtkUtils.py1` look like
editor or manual backup files.

Recommended changes:

- Compare each `.py1` file with its corresponding `.py` file.
- Delete the backup files if they are obsolete.
- If they contain useful changes, merge those changes into the canonical
  source files.
- Add ignore patterns for editor backup files if needed.

### 17. Clean up `.gitignore`

`.gitignore` includes `.github`, even though GitHub workflows are tracked and
should remain visible to contributors.

Recommended changes:

- Remove `.github` from `.gitignore`.
- Keep local-only tool directories such as `.codex`, `.agents`, and `.vscode`
  ignored if they should not be shared.
- Add patterns for generated PyFEM outputs that are not meant to be committed.
- Review ignored files with `git status --ignored` before changing the ignore
  policy.

## Suggested Implementation Order

1. Fix the API naming bug and add a regression test.
2. Correct README links, license metadata, version metadata, and the quick-start
   example.
3. Add `pyproject.toml` with optional dependency groups.
4. Update CI to install the package and run the same test command documented
   for contributors.
5. Fix pytest discovery or standardize fully on unittest.
6. Repair Read the Docs/docs dependency configuration.
7. Replace CLI parsing with `argparse` and test documented options.
8. Remove or quarantine generated outputs and backup files.
9. Tighten documentation builds and link checks after the obvious broken links
   are resolved.

## Verification Commands Used During Audit

```bash
python -m unittest discover -s test -p "*.py"
python -m pytest test
python -m pyfem.core.cli --help
python -m pyfem.core.cli examples/ch02/PatchTest8.pro
```

Observed results:

- `unittest` discovered and ran 189 tests successfully.
- `pytest` collected zero tests.
- `pyfem --help` worked.
- `examples/ch02/PatchTest8.pro` ran successfully.
