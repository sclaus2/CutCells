# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Check that every place that declares the CutCells version agrees.

The version is written in the root pyproject.toml (PyPI), python/pyproject.toml
(builds against an installed core), both CMake projects, the find_package call
of the Python extension, and cutcells.__version__. On a release tag (the
GITHUB_REF environment variable is refs/tags/vX.Y.Z) the tag must match too.
"""

import os
import re
import sys
import tomllib
from pathlib import Path

root = Path(__file__).resolve().parents[2]


def first_match(path, pattern):
    text = (root / path).read_text()
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        sys.exit(f"{path}: no version found (pattern {pattern!r})")
    return match.group(1)


versions = {
    "pyproject.toml": tomllib.loads((root / "pyproject.toml").read_text())["project"]["version"],
    "python/pyproject.toml": tomllib.loads((root / "python/pyproject.toml").read_text())["project"]["version"],
    "python/CMakeLists.txt (project)": first_match(
        "python/CMakeLists.txt", r'^project\(\S+ VERSION "([^"]+)"'
    ),
    "python/CMakeLists.txt (find_package)": first_match(
        "python/CMakeLists.txt", r"find_package\(CutCells ([0-9.]+)"
    ),
    "cpp/CMakeLists.txt": first_match("cpp/CMakeLists.txt", r'^project\(\S+ VERSION "([^"]+)"'),
    "python/cutcells/__init__.py": first_match(
        "python/cutcells/__init__.py", r'^__version__ = "([^"]+)"'
    ),
}

ref = os.environ.get("GITHUB_REF", "")
if ref.startswith("refs/tags/v"):
    versions["git tag"] = ref.removeprefix("refs/tags/v")

for name, version in versions.items():
    print(f"{version:>10}  {name}")

if len(set(versions.values())) != 1:
    sys.exit("The versions above differ: make them agree before releasing.")
print("All versions agree.")
