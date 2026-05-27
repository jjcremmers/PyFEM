"""Sphinx configuration for PyFEM."""

from importlib.metadata import PackageNotFoundError, version as package_version

project = "PyFEM"
copyright = "2026, Joris Remmers"

try:
    release = package_version("pyfem")
except PackageNotFoundError:
    release = "dev"
version = release

root_doc = "index"

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

extensions = [
    "myst_parser",
    "autoapi.extension",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

exclude_patterns = [
    "_build",
    "_api",
    "README.md",
    "img/README.md",
    "Thumbs.db",
    ".DS_Store",
]

autoapi_dirs = ["../pyfem"]
autoapi_root = "api"
autoapi_add_toctree_entry = True
autoapi_keep_files = False
autoapi_ignore = ["*/gui/*"]

autodoc_mock_imports = ["PySide6", "vtk"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

myst_heading_anchors = 3

html_theme = "furo"
html_static_path = ["_static"]
html_favicon = "pyfem.ico"
html_theme_options = {
    "light_logo": "pyfem_logo_official180.png",
    "dark_logo": "pyfem_logo_official180.png",
    "source_repository": "https://github.com/jjcremmers/PyFEM",
    "source_branch": "main",
    "source_directory": "doc/",
}

latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
    "preamble": r"""
        \usepackage{bookmark}
        \setcounter{secnumdepth}{2}
        \setcounter{tocdepth}{3}
    """,
    "figure_align": "htbp",
}
latex_toplevel_sectioning = "chapter"

latex_documents = [
    ("index", "PyFEM.tex", "PyFEM Documentation", "Joris Remmers", "manual"),
]
