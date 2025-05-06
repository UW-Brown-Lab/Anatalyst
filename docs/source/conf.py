# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from unittest.mock import MagicMock

# Add the project root directory to the path so Sphinx can find your modules
sys.path.insert(0, os.path.abspath('../..'))

# Mock modules that might not be installed during documentation generation
# This is helpful for modules with C dependencies or other complex requirements
class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()

MOCK_MODULES = [
    'scanpy', 'scanpy.pp', 'scanpy.tl', 'scanpy.pl', 'scanpy.experimental',
    'numpy', 'pandas', 'matplotlib', 'matplotlib.pyplot', 
    'scipy', 'scipy.sparse', 'scipy.io', 'scipy.stats',
]
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

# -- Project information -----------------------------------------------------
project = 'Brown Lab Downstream Pipeline'
copyright = '2025, Wesley Blashka'
author = 'Wesley Blashka'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',   # Auto-generate documentation from docstrings
    'sphinx.ext.viewcode',  # Add links to source code
    'sphinx.ext.napoleon',  # Support for NumPy and Google style docstrings
    'sphinx.ext.intersphinx',  # Link to other projects' documentation
    'sphinx_rtd_theme',     # Read the Docs theme
    'myst_parser',          # Support for Markdown files
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'

# -- Options for intersphinx extension ---------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable/', None),
}

# -- MyST configuration -----------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
myst_heading_anchors = 3