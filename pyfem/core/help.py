# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

def print_help():
  """
  Print usage instructions for the PyFEM command-line interface.

  This function displays information about how to run PyFEM from the command line,
  available options, a brief description, and a link to the official documentation.
  """

  from pyfem.core.cli import parse_arguments

  try:
    parse_arguments(["--help"])
  except SystemExit:
    pass
