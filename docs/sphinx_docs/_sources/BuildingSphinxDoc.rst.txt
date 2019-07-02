Building Sphinx Doc
===================

Prerequisities
--------------
This documentation uses Sphinx third-party Python package. To generate it from the **docs** folder, you need to install additional components.

Install Python
~~~~~~~~~~~~~~
To install python, you can install for example the `Anaconda package`_ or the `python offical package`_.

Install Sphinx
~~~~~~~~~~~~~~
From the prompt, it is possible to install sphinx, as well as additional libraries required for building this documentation using pip:

.. code-block:: console

	pip install sphinx
	pip install sphinx_rtd_theme
	pip install sphinx-numfig
	pip install sphinxcontrib-bibtex

Building html doc
-----------------
To build the docuentation as a webpage (as shown at: //FP//TODO//inserer adresse), move to the *docs* folder and run under prompt:

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