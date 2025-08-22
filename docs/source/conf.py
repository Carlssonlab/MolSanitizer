# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'MolSanitizer'
copyright = "2024, Carlsson's Lab, Uppsala University"
author = 'Phong Lam, Szymon Pach, Israel Cabeza de Vaca, Ruth Brenk, Jens Carlsson'

release = '0.4.0'
version = '0.4.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    #'sphinx.ext.imgmath',
    'sphinx_copybutton',
    'sphinx.ext.mathjax',
    'sphinx_togglebutton'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']


templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

html_theme_options = {    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'includehidden': True,
    'titles_only': False}

html_static_path = ['_static']
html_css_files = ['custom.css',]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# -- Options for EPUB output
epub_show_urls = 'footnote'
