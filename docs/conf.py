"""Sphinx configuration for the CutCells documentation."""

from pathlib import Path
import sys

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / "python"))

try:
    with (project_root / "python" / "pyproject.toml").open("rb") as f:
        release = tomllib.load(f)["project"]["version"]
except Exception:
    release = "0.4.0"

version = ".".join(release.split(".")[:2])

project = "CutCells"
copyright = "2026, ONERA"
author = "Susanne Claus"

extensions = [
    "sphinx.ext.githubpages",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "sphinx_design",
    "myst_parser",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "__pycache__",
    "quadrays",  # developer notes and plans, not part of the site
]

templates_path = ["_templates"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "logo": {
        "text": "CutCells",
        "alt_text": "CutCells geometry and quadrature documentation",
    },
    "navbar_start": ["navbar-logo"],
    "navbar_center": [],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_persistent": ["search-button"],
    "primary_sidebar_end": [],
    "secondary_sidebar_items": [],
    "show_nav_level": 2,
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sidebar_includehidden": True,
    "show_prev_next": False,
    "search_bar_text": "Search the docs...",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/sclaus2/cutcells",
            "icon": "fab fa-github-square",
            "type": "fontawesome",
        },
    ],
}

html_sidebars = {
    "**": ["sidebar-collapse.html", "global-sidebar.html"],
}

html_title = f"{project} {version}"
htmlhelp_basename = "CutCellsdoc"

pygments_style = "sphinx"
pygments_dark_style = "monokai"

latex_documents = [
    (
        master_doc,
        "CutCells.tex",
        "CutCells Documentation",
        "CutCells Developers",
        "manual",
    ),
]
man_pages = [(master_doc, "cutcells", "CutCells Documentation", [author], 1)]
texinfo_documents = [
    (
        master_doc,
        "CutCells",
        "CutCells Documentation",
        author,
        "CutCells",
        "Cut-cell geometry and quadrature tools",
        "Miscellaneous",
    ),
]
