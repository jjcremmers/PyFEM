# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

"""Unit tests for the public PyFEM API."""

from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory
import unittest

from importlib.metadata import version
import tomllib

import pyfem
from pyfem import run
from pyfem.core.api import PyFEMAPI
from pyfem.core.cli import main as cli_main
from pyfem.core.cli import parse_arguments


class TestPyFEMAPI(unittest.TestCase):
    """Regression tests for the high-level API wrapper."""

    def setUp(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        self.example_dir = repo_root / "examples" / "ch02"
        self.pro_name = "PatchTest8.pro"
        self.dat_name = "PatchTest8.dat"

    def _copy_example(self, target: Path) -> Path:
        copyfile(self.example_dir / self.pro_name, target / self.pro_name)
        copyfile(self.example_dir / self.dat_name, target / self.dat_name)
        return target / self.pro_name

    def test_run_returns_results_dict(self) -> None:
        """Test that the top-level run helper completes an analysis."""
        with TemporaryDirectory() as tmp:
            pro_file = self._copy_example(Path(tmp))

            results = run(pro_file)

            self.assertFalse(results["active"])
            self.assertIn("globdat", results)
            self.assertIn("props", results)

    def test_snake_case_and_camel_case_api_names_work(self) -> None:
        """Test current snake_case API and backward-compatible camelCase names."""
        with TemporaryDirectory() as tmp:
            pro_file = self._copy_example(Path(tmp))
            api = PyFEMAPI(pro_file)

            self.assertTrue(api.is_active)
            self.assertTrue(api.isActive)

            api.run_all()
            self.assertFalse(api.is_active)

            results = api.get_results()
            legacy_results = api.getResults()
            self.assertIs(results["globdat"], legacy_results["globdat"])

            api.runAll()
            api.close()


class TestPackageMetadata(unittest.TestCase):
    """Tests for public package metadata."""

    def test_package_version_matches_installed_metadata(self) -> None:
        """Test that pyfem.__version__ follows the installed package version."""
        self.assertEqual(pyfem.__version__, "2026.9")
        self.assertEqual(version("pyfem"), pyfem.__version__)

    def test_gui_dependency_is_optional_and_vtk_is_core(self) -> None:
        """Test that VTK is installed by default while GUI remains optional."""
        repo_root = Path(__file__).resolve().parents[1]
        pyproject = tomllib.loads((repo_root / "pyproject.toml").read_text())

        dependencies = set(pyproject["project"]["dependencies"])
        optional = pyproject["project"]["optional-dependencies"]

        self.assertNotIn("PySide6", dependencies)
        self.assertIn("vtk", dependencies)
        self.assertEqual(optional["gui"], ["PySide6"])
        self.assertIn("myst-parser>=2.0.0", optional["docs"])
        self.assertIn("sphinx-rtd-theme>=1.2.0", optional["docs"])
        self.assertNotIn("vtk", optional["docs"])
        self.assertIn("pytest", optional["dev"])
        self.assertIn("coverage", optional["dev"])
        self.assertIn("build", optional["dev"])
        self.assertEqual(optional["all"], ["PySide6"])


class TestCommandLineInterface(unittest.TestCase):
    """Tests for command-line argument parsing and execution."""

    def setUp(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        self.example_dir = repo_root / "examples" / "ch02"
        self.pro_name = "PatchTest8.pro"
        self.dat_name = "PatchTest8.dat"

    def _copy_example(self, target: Path) -> Path:
        copyfile(self.example_dir / self.pro_name, target / self.pro_name)
        copyfile(self.example_dir / self.dat_name, target / self.dat_name)
        return target / self.pro_name

    def test_parse_positional_input(self) -> None:
        """Test positional input-file parsing."""
        args = parse_arguments(["model.pro"])

        self.assertEqual(args.input_file, "model.pro")
        self.assertIsNone(args.dump_file)
        self.assertEqual(args.parameters, [])

    def test_parse_input_option_and_parameter_overrides(self) -> None:
        """Test -i, -d, and repeated -p options."""
        args = parse_arguments(
            ["-i", "model.pro", "-d", "state.dump", "-p", "E=210000", "-p", "nu=0.3"]
        )

        self.assertEqual(args.input_file, "model.pro")
        self.assertEqual(args.dump_file, "state.dump")
        self.assertEqual(args.parameters, ["E=210000", "nu=0.3"])

    def test_parse_rejects_duplicate_input_forms(self) -> None:
        """Test that positional and -i input forms are mutually exclusive."""
        with self.assertRaises(SystemExit):
            parse_arguments(["model.pro", "-i", "other.pro"])

    def test_parse_version_exits_successfully(self) -> None:
        """Test --version handling."""
        with self.assertRaises(SystemExit) as cm:
            parse_arguments(["--version"])

        self.assertEqual(cm.exception.code, 0)

    def test_cli_runs_input_option(self) -> None:
        """Test that the documented -i input option runs an analysis."""
        with TemporaryDirectory() as tmp:
            pro_file = self._copy_example(Path(tmp))

            cli_main(["-i", str(pro_file)])


if __name__ == "__main__":
    unittest.main()
