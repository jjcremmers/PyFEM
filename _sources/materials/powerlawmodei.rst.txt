================
PowerLawModeI
================

``PowerLawModeI`` defines a mode-I cohesive traction–separation relation using
fracture energy and peak traction.

Overview
--------

Material type: ``PowerLawModeI``

Parameters
----------

Mandatory Parameters
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Must be set to ``"PowerLawModeI"``
   * - ``Gc``
     - Fracture energy
   * - ``Tult``
     - Ultimate traction

Examples
--------

- Interface modeling in cohesive elements

See Also
--------
- :doc:`../elements/interface`
