Compiling Anima from source
===========================

Requirements
^^^^^^^^^^^^
Anima is multi-platform and will compile and run on the major platforms: Windows, OSX and Linux. It requires two major things: 

* `CMake <http://www.cmake.org>`_ (>= 3.1.0) as a cross-platform Makefile generation software (we use its super-project ability to also download and compile its dependencies)
* A compilation environment: Visual studio 2015 on Windows, Xcode and developer tools on OSX, gcc/g++ or any other C++ compiler on Linux

Code architecture
^^^^^^^^^^^^^^^^^

Anima itself follows a modular structure, organized in 6 main modules, each one associated to a folder on the main repository:

* Math-tools. This module contains generic mathematical tools: statistical tests, statistical distributions, optimizers, spherical harmonics...
* Filtering. This module contains image filters: smoothing, denoising...
* Diffusion. This module includes diffusion model estimation as well as tractography algorithms.
* Registration. This one holds registration tools: resamplers, interpolators, transformations, registration algorithms...
* Segmentation. segmentation tools
* Quantitative MRI. Estimation of quantitative parameters from MRI, MR simulation

Some of these modules are dependent on others (e.g. math-tools is the base to everything). If you know what you are doing, you can choose to compile only a subset of these.

Compilation Instructions
^^^^^^^^^^^^^^^^^^^^^^^^

Anima is available as a superproject, including all its modules and links and instructions for its dependencies. Building ANIMA is simple:

* clone the repository from github (use the first line by default, the second if you have set up your SSH keys): 

.. code-block:: sh

	git clone https://github.com/Inria-Visages/Anima-Public.git src
	git clone git@github.com:Inria-Visages/Anima-Public.git src

* then, run CMake in a new build folder, change any options if you wish to change the default compilation (which downloads and compiles all dependencies and tools): 

.. code-block:: sh

	mkdir build
	cd build
	ccmake ../src

* build using your environment (a `make` or `ninja` will be enough on Linux and OSX, open Visual Studio on Windows)

And there you go, the `bin` folder in the build directory will contain all tools described in this documentation.
