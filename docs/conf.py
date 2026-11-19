import pathlib
import re

project = "StdFace Documentation"
author = "ISSP Center"

# The C generator's version lives in src/version.h only; read it from there.
_header = (pathlib.Path(__file__).resolve().parent.parent / "src" / "version.h").read_text()
_number = dict(re.findall(r'#define\s+VERSION_(MAJOR|MINOR|PATCH)\s+(\d+)', _header))
_prerelease = re.search(r'#define\s+VERSION_PRERELEASE\s+"([^"]*)"', _header).group(1)
version = "{MAJOR}.{MINOR}".format(**_number)
release = "{MAJOR}.{MINOR}.{PATCH}".format(**_number)
if _prerelease:
    release += "-" + _prerelease

# The Python package is versioned independently, in its version module only.
_py_module = (pathlib.Path(__file__).resolve().parent.parent
              / "python" / "stdface" / "core" / "version.py").read_text()
_py_version = re.search(r'^__version__\s*=\s*"([^"]+)"', _py_module, re.M).group(1)
rst_epilog = f".. |py_version| replace:: {_py_version}\n"
extensions = []
templates_path = ['_templates']
exclude_patterns = ['_meta', '._*', '_build']
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']
master_doc = 'index'
