Installing Anima scripts
========================

Requirements
------------

Anima scripts does not require many things. But here is a list of what is required:

* `Python <https://www.python.org>`_ (>= 3.0): for most scripts, Python is required. Python 2 or 3 are supported although python 3 and more recent is advised
* numpy, scipy, pandas and pydicom are required for some scripts. The following command should do the trick to install them:

.. code-block:: sh

	pip3 install numpy pydicom pandas scipy --user

* Anima: either built from source or the releases. Some scripts (e.g. in diffusion) require non public code that are in the releases but not the open-source code unless you have access to our repositories. See there how to :doc:`compile it from source <compile_source>` or :doc:`install the binaries <install_binaries>`
* A cluster with an `OAR <http://oar.imag.fr>`_ scheduler for the atlasing scripts

Anima scripts installation instructions
---------------------------------------

* Install Git LFS if you want to clone the additional data repository (if not, skip this step). It may be downloaded from `here <https://git-lfs.github.com/>`_. Follow then the installation instructions (on Mac and linux, simply run the install.sh script as a super-user)

* Get the additional data repository by cloning it or from our `website <https://anima.irisa.fr/downloads/>`_.

.. code-block:: sh

	git clone https://github.com/Inria-Visages/Anima-Scripts-Data-Public.git

* Get the release of Anima scripts from our `website <https://anima.irisa.fr/downloads/>`_ or clone the repository using 

.. code-block:: sh

	git clone https://github.com/Inria-Visages/Anima-Scripts-Public.git

* Run the configure.py file inside Anima scripts public folder with your configuration. The default values assume you have followed this guide. An example of use:

.. code-block:: sh

	python3 configure.py -a ~/Anima-Public/build/bin -s ~/Anima-Scripts-Public -d ~/Anima-Scripts-Data-Public

* The options are as follows:
  * ``-s`` should link to where you cloned the Anima scripts repository
  * ``-a`` should point to where you have the Anima executables
  * ``-d`` should point to where you have cloned the additional data repository

And now you are all setup to go, you can now read about all scripts in the other sections of the documentation.
