# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'MolSanitizer'
copyright = "2026, Carlsson's Lab, Uppsala University"
author = 'Phong Lam, Szymon Pach, Israel Cabeza de Vaca, Ruth Brenk, Jens Carlsson'

release = "0.7.3rc1"
version = "0.7.3rc1"

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # For Google/NumPy docstrings
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_copybutton',
    'sphinx.ext.mathjax',
    'sphinx_togglebutton'
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']


templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_logo = '../../assets/icon_text.png'

html_theme_options = {    # Toc options
    'logo_only': True,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'includehidden': True,
    'titles_only': False}

html_static_path = ['_static']
html_css_files = ['custom.css',]

# The suffix(es) of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# -- Options for EPUB output
epub_show_urls = 'footnote'
