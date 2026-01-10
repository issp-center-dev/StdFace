API Reference
=============

This section documents the public API defined in headers under ``src/``.

.. note::

   Unlike typical C projects, StdFace does not have an ``include/`` directory.
   All public headers are located in ``src/``.

Headers
-------

.. toctree::
   :maxdepth: 1

   StdFace_ModelUtil
   StdFace_vals
   setmemory
   export_wannier90
   version

Summary
-------

+------------------------+--------------------------------------------------+
| Header                 | Description                                      |
+========================+==================================================+
| ``StdFace_ModelUtil.h``| Core model utilities (37 functions)              |
+------------------------+--------------------------------------------------+
| ``StdFace_vals.h``     | Configuration structure StdIntList (~125 fields) |
+------------------------+--------------------------------------------------+
| ``setmemory.h``        | Memory allocation utilities (24 functions)       |
+------------------------+--------------------------------------------------+
| ``export_wannier90.h`` | Wannier90 export (HWAVE only, 2 functions)       |
+------------------------+--------------------------------------------------+
| ``version.h``          | Version information (1 function, 4 macros)       |
+------------------------+--------------------------------------------------+
