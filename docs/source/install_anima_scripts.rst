Installing Anima scripts
========================

Requirements
------------

Anima scripts does not require many things. But here is a list of what is required:

* `Python <https://www.python.org>`_ (>= 3.0): for most scripts, Python is required. Python 2 or 3 are supported although python 3 and more recent is advised
* numpy, scipy, pandas and pydicom are required for some scripts. The following command should do the trick to install them:

.. prompt:: bash $

	pip install numpy pydicom pandas scipy

* Anima: either built from source or the releases. Some scripts (e.g. in diffusion) require non public code that are in the releases but not the open-source code unless you have access to our repositories. See there how to :doc:`compile it from source <compile_source>` or :doc:`install the binaries <install_binaries>`
* `Anima scripts data <https://github.com/Inria-Visages/Anima-Scripts-Data-Public>`_: data required for some scripts to work (brain extraction and diffusion scripts)
* A cluster with an `OAR <http://oar.imag.fr>`_ scheduler for the atlasing scripts

Anima scripts installation instructions
---------------------------------------

* Get the additional data repository by cloning it. Be careful, you should have installed `Git LFS <https://git-lfs.github.com/>`_ first to be able to get it

.. prompt:: bash $

	git clone https://github.com/Inria-Visages/Anima-Scripts-Data-Public.git

* Get the release of Anima scripts from our `website <https://inria-visages.github.io/Anima-Public/downloads>`_ or clone the repository using 

.. prompt:: bash $

	git clone https://github.com/Inria-Visages/Anima-Scripts-Public.git

* Copy / paste the **example-config.txt** file in **.anima/config.txt** in your home folder
* Update the paths in your config file to match where Anima binaries are and where the additional data folder is (use full paths, no tilde)

  * **anima-scripts-public-root** should link to where you cloned this repository
  * **anima** should point to where you have the Anima executables
  * **extra-data-root** should point to where you have cloned the additional data repository

And now you are all setup to go, you can now read about all scripts in the other sections of the documentation.
