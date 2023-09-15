# -- Dockerfile --
# This Dockerfile provides a simple way to build documentation
# without requiring to install all dependencies locally

# Pull Ubuntu latest version
FROM ubuntu:20.04

# Minimal Python installation
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y python3 \
                    python3-pip

# Install Sphinx
RUN pip3 install sphinx sphinx_rtd_theme sphinx-numfig sphinxcontrib-bibtex

# # Install LaTeX (to build PDF version of the doc)
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y texlive \
                    texlive-latex-extra \
                    latexmk

WORKDIR /home
ENTRYPOINT ["make"]
CMD ["html"] 

# -- User guide --

# To build the container image using this Dockerfile
# docker image build -t ecogen/sphinx-doc .

# To run the container 
# For HTML output
# docker container run --mount type=bind,source=${PWD},target=/home ecogen/sphinx-doc
# Or 
# docker container run --mount type=bind,source=${PWD},target=/home ecogen/sphinx-doc html
# For PDF output
# docker container run --mount type=bind,source=${PWD},target=/home ecogen/sphinx-doc pdf

# To troubleshoot any bug it is also possible to run the container interactively
# docker container run --mount type=bind,source=${PWD},target=/home -it -w /home --entrypoint /bin/bash ecogen/sphinx-doc
# make html