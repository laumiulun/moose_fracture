# MOOSE Flavored Markdown

MOOSE generates documentation using the [python-markdown] package,
which includes the ability to use extensions from others as well as define custom extensions. The following pages outline the
outside extensions enabled, as well as the custom extensions included with [MOOSE],
including the custom syntax defined exclusively for documenting [MOOSE] source code.

## MooseDocs Extensions
Much of the capability that enables the MooseDocs system to be useful for creating documentation
stems from the set of custom extensions created. In general, these extensions add commands, which
begin with an exclamation point (!) and some number of arguments. Therefore, each extension contains
configuration options at the extension level and each command defined as a settings that can be
applied.

\ref{moose-extensions} summarizes all of the available MooseDocs extensions and provides a link to
another page that details the use and configuration of the extension as well as a list of any
commands that add associated settings defined in the extension.

!table id=moose-extensions caption=List of MooseDocs extensions for writing documentation.
| Name | Description |
| - | - |
| [global](extensions/global.md) | Implements a global list of Markdown hyperlinks. |
| [bibtex](extensions/bibtex.md) | Allows LaTeX/BibTeX style citations and bibliographies. |
| [css](extensions/css.md) | Enables [CSS](https://en.wikipedia.org/wiki/Cascading_Style_Sheets) settings to be applied via markdown |
| [devel](extensions/devel.md) | Adds the ability to extract configuration and commands settings in addition to build status and MOOSE package information (MOOSE developer tools). |
| [include](extensions/include.md) | Adds ability to include complete or partial markdown files. |
| [misc](extensions/misc.md) | Adds code copy button and scrolling contents to html. |
| [media](extensions/media.md) | Adds markdown for including images, sliders, and videos. |
| [tables](extensions/tables.md) | Markdown syntax for numbered tables. |
| [listings](extensions/listings.md) | Markdown syntax for numbered code blocks. |
| [refs](extensions/refs.md) | Implements latex style references to numbered floats (e.g., figures and tables) |
| [templates](extensions/templates.md) | Allows for converted markdown to be applied to arbitrary templates. |
| [app_syntax](extensions/app_syntax.md) | Adds markdown syntax for extracting content from MOOSE source code. |
| [sqa](extensions/sqa.md) | Adds tools for creating software quality documentation using MOOSE templates. |
| [gchart](extensions/gchart.md) | Adds [Google chart](https://developers.google.com/chart/) capability to MOOSE markdown. |
| [katex](extensions/katex.md) | Adds math support via [KaTeX](https://khan.github.io/KaTeX/). |


## Python-Markdown Extensions

The [python-markdown] package includes many useful, officially-supported extensions, as listed on the
[available extensions](https://pythonhosted.org/Markdown/extensions/) page.
\ref{official-extensions} is a list the extensions that are utilized on this website and will
likely be of use as you develop your own  website, reports, or presentations.

!table id=official-extensions caption=List of official [python-markdown] extensions useful for writing documentation.
| Name | Description |
| - | - |
| [toc](https://pythonhosted.org/Markdown/extensions/toc.html) | Generates a Table of Contents from a Markdown document and adds it into the resulting HTML document. |
| [smarty](https://pythonhosted.org/Markdown/extensions/smarty.html) | The SmartyPants extension converts ASCII dashes, quotes and ellipses to their HTML entity equivalents. |
| [admonition](https://pythonhosted.org/Markdown/extensions/admonition.html) | The Admonition extension adds [rST-style](http://docutils.sourceforge.net/docs/ref/rst/directives.html#specific-admonitionss) admonitions to Markdown documents. |
| [extra](https://pythonhosted.org/Markdown/extensions/extra.html) | A compilation of various Python-Markdown extensions that (mostly) imitates [PHP Markdown Extra](https://michelf.ca/projects/php-markdown/extra/). |
| [meta](https://pythonhosted.org/Markdown/extensions/meta_data.html) | Adds a syntax for defining meta-data about a document, which is used by MooseDocs to specify template arguments in the [MooseDocs.extensions.template](extensions/templates.md) extension. |
| [mdx_math](https://github.com/mitya57/python-markdown-math) | This extension adds math formulas support to Python-Markdown. |

The [admonition](https://pythonhosted.org/Markdown/extensions/admonition.html) package enables for important and critical
items to be highlighted, using the syntax detailed below and the package documentation: [admonition](https://pythonhosted.org/Markdown/extensions/admonition.html).

```markdown
!!! type "An optional title"
    A detailed message paragraph that is indented by 4 spaces and can include any number of lines.
```

The supported "types" for MOOSE are: "info", "note", "important, "warning", "danger", and "error":

!!! info "Optional Info Title"
    This is some information you want people to know about.

!!! note "Optional Note Title"
    This is an example of a note.

!!! important "Optional Important Title"
    This is an example of something important.

!!! warning "Optional Warning Title"
    This is a warning.

!!! danger "Optional Danger Title"
    This is something very dangerous.

!!! error "Optional Error Title"
    This is an error message.
