=====
Dummy
=====

``Dummy`` is a simple linear traction law for interface testing. It uses a
scalar stiffness `D` and supports 2D or 3D traction outputs.

Overview
--------

Material type: ``Dummy``

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
     - Must be set to ``"Dummy"``
   * - ``D``
     - Interface stiffness

Examples
--------

- Simple interface tests and validation
