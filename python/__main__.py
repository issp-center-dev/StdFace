"""Entry point for ``python python/__main__.py``.

Delegates to ``stdface.__main__.main()``.
"""
from stdface.__main__ import main
import sys

if __name__ == "__main__":
    sys.exit(main())
