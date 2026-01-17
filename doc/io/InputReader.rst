===========
InputReader
===========

``InputReader`` reads a ``.pro`` input file and constructs the analysis state
(``props`` and ``globdat``). It can also read a previously saved dump to
restore the full state.

Overview
--------

Function: ``InputReader(argv)`` and ``InputRead(fname, dname=None, parameters=None)``

- Reads properties from a ``.pro`` file; parses mesh, elements, dofs, models.
- Sets up logging, prefix, and initial global state.
- Optionally reads a dump file (pickle) to restore state instead of parsing.

Command-Line Arguments
----------------------

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Option
     - Description
   * - ``-i``, ``--input``
     - Path to the ``.pro`` file
   * - ``-d``, ``--dump``
     - Path to a pickle dump (created by ``DataDump``) to restore state
   * - ``-p``, ``--param``
     - Override parameters as ``name=value`` (can be repeated)
   * - ``-h``, ``--help``
     - Show help

Programmatic Use
----------------

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Function
     - Description
   * - ``InputRead(fname, dname=None, parameters=None)``
     - Returns ``props, globdat``; reads from ``fname`` or ``dname`` when provided.

Examples
--------

- Read from a ``.pro`` file:

.. code-block:: python

   from pyfem.io.InputReader import InputRead
   props, globdat = InputRead("examples/ch03/cantilever8.pro")

- Restore from a dump:

.. code-block:: python

   props, globdat = InputRead(None, dname="results/run.dump")

See Also
--------
- :doc:`overview`
- Related: :doc:`DataDump`
