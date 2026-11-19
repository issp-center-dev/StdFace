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
extensions = []
templates_path = ['_templates']
exclude_patterns = ['_meta', '._*', '_build']
html_theme = 'alabaster'
html_static_path = []
master_doc = 'index'
