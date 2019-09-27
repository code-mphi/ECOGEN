Building Sphinx Doc
===================

Prerequisities
--------------
This documentation uses Sphinx third-party Python package. To generate it locally from the **docs/sphinx_docs/** folder, you need to install additional components.

Install Python
~~~~~~~~~~~~~~
To install python, you can install for example the `Anaconda package`_ or the `python offical package`_.
You can also install the packages directly from ubuntu terminal:

.. code-block:: console

	apt-get install python3
	apt-get install python-pip

Install Sphinx
~~~~~~~~~~~~~~
From the prompt, it is possible to install sphinx, as well as additional libraries required for building this documentation using pip:

.. code-block:: console

	pip install sphinx
	pip install sphinx_rtd_theme
	pip install sphinx-numfig
	pip install sphinxcontrib-bibtex

Install Latex
~~~~~~~~~~~~~
If you want to build pdf, you will need Latex installed

.. code-block:: console

	apt-get install texlive
	apt-get install texlive-latex-extra
	apt-get install latexmk

Building html doc
-----------------
To build the docuentation as a webpage (as shown on ECOGEN_ webSite), move to the *docs* folder and run under prompt:

.. code-block:: console

	make html

Building pdf doc
----------------
To build the docuentation as a PDF, move to the *docs* folder and run under prompt:

.. code-block:: console

	make latexpdf

This will generate a folder *docs/build/latex* containing source files in *LaTeX* that can be used to generate a pdf.

Learning Sphinx
---------------
To learn how to developp a documentation using sphinx, here are some usefull links:

- `Sphinx documentation`_
- Hosting documentation and read the docs theme : `Read the docs website`_


.. _`Anaconda package`: https://www.anaconda.com/distribution/
.. _`python offical package` : https://www.python.org/
.. _`Sphinx documentation`: https://www.sphinx-doc.org/en/master/contents.html
.. _`Read the docs website`: https://readthedocs.org/
.. _ECOGEN: https://code-mphi.github.io/ECOGEN/docs/sphinx_docs/index.html