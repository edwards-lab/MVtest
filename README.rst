MVtest GWAS Analysis
====================

MVtest is an analysis tool that can be run on many common file formats using
syntax similar to programs you've probably already used. Instructions for
installation can be found below.

MVtest dependencies are currently incompatible with Python 3.12
===============================================================

As of May 2024, we have been made aware that a dependency of a library we rely on
for parsing bgen data has a dependency of its own which uses a deprecated feature
that makes it incompatible with Python 3.12 (and greater). Until this has been 
fixed, MVtest can only be installed on systems that run on Python 3.10 (or earlier
until 3.7). 

Installation
============

MVtest requires python 3.7.x as well as the following libraries:

* NumPy (version 1.16.2 or later)           www.numpy.org
* SciPY (version 1.3.0 or later)            www.scipy.org
* pytabix (version 0.1 or later)            https://pypi.org/project/pytabix/
* bgen-reader (version 3.0.6 or later)      https://pypi.org/project/bgen-reader/
* libGWAS (version 1.2.0 or later)          https://github.com/edwards-lab/libGWAS

MVtest's installation will attempt to install these required components
for you, however, it requires that you have write permission to the
installation directory. If you are using a shared system and lack the
necessary privileges to install libraries and software yourself, you should
please see one of the sections, Miniconda_ or virtual-env_ below
for instructions on different options for setting up your own python
environement which will exist entirely under your own control.

Installation
===================

To install MVTest, simply clone the sources using the following command:


$ `git clone https://github.com/edwards-lab/MVtest`

Or you may visit the website and download the tarball directly from github: https://github.com/edwards-lab/MVtest

Once you have downloaded the software, simply extract the contents and run the
following command to install it:

$ `pip install .`

If no errors are reported, it should be installed and ready to use.

**Regarding PYTHON 2** I have completely switched over to python 3 without 
trying to remain compatible with Python 2 because the bgen_reader no longer
supports 2 and the end of life is only a few months from my writing this. 

As such, if you wish to use python2, you will need to install an older version.

System Requirements
+++++++++++++++++++
Aside from the library dependencies, MVtest's requirements depend largely on
the number of SNPs and individuals being analyzed as well as the data format
being used. In general, GWAS sized datasets will require several gigabytes of
memory when using the traditional pedigree format, however, even 10s of
thousands of subjects can be analyzed with less than 1 gigabyte of RAM when
the data is formatted as transposed pedigree or PLINK's default bed format.

Otherwise, it is recommended that the system be run on a unix-like system
such as Linux or OS X, but it should work under windows as well (we can't
offer support for running MVtest under windows).

Running Unit Tests
++++++++++++++++++
MVtest comes with a unit test suite which can be run prior to installation.
To run the tests, simply run the following command from within the root
directory of the extracted archive's contents:

$ `pytest`

If no errors are reported, then mvtest should run correctly on your system.

.. _virtual-env:

Virtual Env
+++++++++++
Virtual ENV is a powerful too for python programmers and end users alike as it
allows for users to deploy different versions of python applications without
the need for root access to the machine.

Because MVtest requires version 3.7, you'll need to ensure that your machine's
python version is in compliance. Virtual Env basically uses the the system
version of python, but creates a user owned environment wrapper allowing
users to install libraries easily without administrative rights to the
machine.

For a helpful introduction to VirtualEnv, please have a look at the
tutorial: http://www.simononsoftware.com/virtualenv-tutorial/

.. _miniconda:

Miniconda
+++++++++
Miniconda is a minimal version of the package manager used by the Anaconda
python distribution. It makes it easy to create local installations of python
with the latest versions of the common scientific libraries for users who don't
have root access to their target machines. Basically, when you use miniconda,
you'll be installing your own version of Python into a directory under your
control which allows you to install anything else you need without having to
submit a helpdesk ticket for administrative assistance.

Unlike pip, the folks behind the conda distributions provide binary downloads
of it's selected library components. As such, only the most popular libraries,
such as pip, NumPY and SciPy, are supported by conda itself. However, these do
not require compilation and may be easier to get installed than when using
pip alone. I have experienced difficulty installing SciPy through pip and setup
tools on our cluster here at vanderbilt due to non-standard paths for certain
required components, but mini-conda always comes through.

Firstly, download and install the appropriate version of miniconda at the
project website. Please be sure to choose the Python 2 version:
http://conda.pydata.org/miniconda.html

While it is doing the installation, please allow it to update your PATH
information. If you prefer not to always use this version of python in the
future, simple tell it not to update your .bashrc file and note the
instructions for loading and unloading your new python environment. Please
note that even if you chose to update your .bashrc file, you will need to
follow directions for loading the changes into your current shell.

Once those changes have taken effect, install setuptools and scipy:
$ `conda install pip scipy`

Installing SciPy will also force the installation of NumPy, which is
also required for running mvtest. (setuptools includes easy_install).

Once that has been completed successfully, you should be ready to follow
the standard instructions for installing mvtest.


MVtest Online Manual
====================

The online manual can be found at http://edwards-lab.github.io/MVtest/

For developers who would like to use the GWAS parsers, the API manual can be
found at http://edwards-lab.github.io/MVtest/api/index.html
