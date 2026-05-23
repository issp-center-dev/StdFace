Input File Format
=================

Basic Syntax
------------

StdFace input files use a simple ``key = value`` format:

.. code-block:: text

   model = "Hubbard"
   lattice = square
   W = 2
   L = 2
   t = 1.0
   U = 1.0

Formatting rules:

- One key-value pair per line
- Keys and values separated by ``=``
- Whitespace around ``=`` is optional (both ``W=2`` and ``W = 2`` are accepted)
- String values may be quoted (``"Hubbard"``) or unquoted (``square``)
- Numeric values are written without quotes
- Keywords and most string values are case-insensitive; they are normalised to
  lowercase internally (e.g., ``model = "Hubbard"`` and ``model = "hubbard"``
  are equivalent)

Comment Syntax
--------------

Lines beginning with ``//`` are treated as comments and ignored:

.. code-block:: text

   // This is a comment line
   model = "Hubbard"
   lattice = square

Empty lines are also skipped.  **Inline comments are not supported** — the
parser splits each line on the first ``=`` sign, so a trailing ``// comment``
would be included in the value and cause a parse error.

Ordering
--------

Key ordering within the input file is not required.  Parameters may appear
in any order.

