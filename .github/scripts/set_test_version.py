# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Give the Python package a version of its own for an upload to TestPyPI.

    python set_test_version.py 0.6.0rc1

TestPyPI, like PyPI, never takes the same version twice, so each trial upload
needs a new one. Only the Python metadata changes (both pyproject.toml files
and cutcells.__version__); the CMake projects keep their numeric version, which
a pre-release such as 0.6.0rc1 could not be.
"""

import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[2]

if len(sys.argv) != 2 or not sys.argv[1]:
    sys.exit("Give the version for TestPyPI, e.g. 0.6.0rc1 or 0.5.1.dev1.")
version = sys.argv[1]
# public PEP 440 versions: release, then optional pre-, post- and dev-release
if not re.fullmatch(r"\d+(\.\d+)*((a|b|rc)\d+)?(\.post\d+)?(\.dev\d+)?", version):
    sys.exit(f"{version!r} is not a PEP 440 version such as 0.6.0rc1 or 0.5.1.dev1.")

for path, pattern in [
    ("pyproject.toml", r'^version = "[^"]+"'),
    ("python/pyproject.toml", r'^version = "[^"]+"'),
    ("python/cutcells/__init__.py", r'^__version__ = "[^"]+"'),
]:
    file = root / path
    text, n = re.subn(pattern, lambda m: m.group(0).split("=")[0] + f'= "{version}"', file.read_text(),
                      count=1, flags=re.MULTILINE)
    if n != 1:
        sys.exit(f"{path}: no version line found")
    file.write_text(text)
    print(f"{version:>12}  {path}")
