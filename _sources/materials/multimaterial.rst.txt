=============
MultiMaterial
=============

The ``MultiMaterial`` model groups multiple material models for layered or
zone-based analyses. It delegates `getStress()` to the active sub-material
specified by the element (e.g., via `iMat`).

Overview
--------

Material type: ``MultiMaterial``

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
     - Must be set to ``"MultiMaterial"``
   * - ``materials``
     - List of material names referenced as child blocks

Optional Parameters
~~~~~~~~~~~~~~~~~~~

Child material blocks must define their own parameters. Each child is any valid
material model (e.g., ``Isotropic``, ``TransverseIsotropic``).

Examples
--------

- Layered laminates with SLS elements

See Also
--------
- :doc:`../elements/sls`
