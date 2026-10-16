# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CutCells'
copyright = '2024, ONERA'
author = 'Susanne Claus'
release = '0.2.0'
import cutcells

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.imgmath',
              'sphinx.ext.autosummary',
              'sphinx.ext.ifconfig',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'sphinx.ext.todo',
              'sphinx_copybutton',
#              "autoapi.extension",
              ]

templates_path = ['_templates']
exclude_patterns = []

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = False

autoclass_content = 'both'  # this concatenates the class.doc + __init__.doc
autodoc_typehints = "description"  # typehints are nicely placed

todo_include_todos = True
todo_link_only = True

# alternatively, could use autoapi:
# # autoapi coonfiguration:
# package_dir = next(x for x in nanobind_example.__path__ if "site-packages" in x)
# autoapi_dirs = [f"{package_dir}/"]
# # https://github.com/readthedocs/sphinx-autoapi/issues/405
# #autoapi_ignore = ["*__init__.py*"]
# autoapi_member_order = "groupwise"
# # autoapi_generate_api_docs = False



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'furo'
html_theme = 'nature'
html_static_path = ['_static']
