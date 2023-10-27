# ihop documentation build configuration file
# Configuration file for the Sphinx documentation builder. Oct 27 2023.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Extension information ---------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('_extensions'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = u'ihop'
copyright = u'2023, ihop contributors'
author = u'Ivana Escobar'


# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#

from subprocess import check_output, CalledProcessError

def get_version():
    """
    Return the latest tag (checkpoint) and, if there have
    been commits since the version was tagged, the commit hash.

    To get just the release tag use:
    version = version.split('-')[0]
    """

    try:
        version = check_output(['git', 'describe', '--tags', '--always'],
                               universal_newlines=True)
    except CalledProcessError:
        return 'unknown version'

    return version.rstrip()

# "version" is used for html build
version = get_version()
# "release" is used for LaTeX build
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinxcontrib.bibtex',
              'sphinxcontrib.programoutput',
              'ihop']

autodoc_mock_imports = ['matplotlib', 'mpl_toolkits']

bibtex_bibfiles = ['ihop_references.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. 
# Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Do not highlight code blocks unless a language is specified explicitly.
highlight_language = 'none'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# MITgcm: number figures

numfig_format = {'figure': 'Figure %s',
                 'table': 'Table %s',
                 'code-block': 'Code %s',
                }
numfig = True

# number figures within section
numfig_secnum_depth = 1

#math_number_all = True
numfig_format = {'figure': 'Figure %s', 'table': 'Table %s', 'code-block': 'Listing %s', 'section': 'Section %s'}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'renku'

html_static_path = ['_static']

html_css_files = [
    'css/custom.css',
    'css/wrap_tables.css',
]

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    'preamble': r'''
    \setcounter{secnumdepth}{3}
    \newcommand{\p}[1]{\frac{\partial }{\partial #1}}
    \newcommand{\pp}[2]{\frac{\partial #1}{\partial #2}}
    \newcommand{\dd}[2]{\frac{d #1}{d #2}}
    \newcommand{\h}{\frac{1}{2}}
    \setlength{\tymax}{0.5\textwidth}
    ''',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'ihop.tex', u'ihop Documentation',
     u'Ivana Escobar', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'ihop', u'ihop Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'ihop', u'ihop Documentation',
     author, 'ihop', 'A package for underwater acoustics for MITgcm.',
     'Miscellaneous'),
]
