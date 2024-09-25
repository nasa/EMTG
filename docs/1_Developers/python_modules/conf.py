# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../../PyEMTG/ConvergenceAnimator/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/Converters/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/GMAT_mission_renderer/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/HighFidelity/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/JourneyOptionsPanel/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/OptionsOverhaul/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/PEATSA/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/PyGMATscripter/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/ReportGenerators/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/SimpleMonteCarlo/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/SpacecraftOptionsPanel/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/SpiceyPy_Utilities/'))
sys.path.insert(0, os.path.abspath('../../../PyEMTG/'))


# -- Project information -----------------------------------------------------

project = 'EMTG python'
copyright = '2024, Jacob Englander, Donald Ellison, Jeremy Knittel, Noble Hatten'
author = 'Jacob Englander, Donald Ellison, Jeremy Knittel, Noble Hatten'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
'sphinx.ext.autodoc',
'sphinx.ext.viewcode',
'sphinx.ext.napoleon'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_autogendocs', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If true, the reST sources are included in the HTML build as _sources/name. The default is True.
html_copy_source = False