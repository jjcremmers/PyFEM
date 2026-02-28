# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

def print_help():
  """
  Print usage instructions for the PyFEM command-line interface.

  This function displays information about how to run PyFEM from the command line,
  available options, a brief description, and a link to the official documentation.
  """

  print("""
PyFEM Command Line Usage:

  pyfem <inputfile>

Options:
  --help, -h        Show this help message and exit

Description:
  Run a PyFEM analysis using the specified input file.

For more information and documentation, visit:
  https://jorisremmers.com/PyFEM
""")
