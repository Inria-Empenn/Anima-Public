Anima documentation
===================

Welcome to the documentation of Anima and Anima scripts. Thanks for your interest. Here are a few directies to installing and what Anima and Anima scripts may do for you. Mailing lists are available for any questions on the tools at anima-users@inria.fr

First steps
-----------

If this is your first time here, you may be very well interested in installing Anima. This may be done either from :doc:`source compilation <compile_source>` or by :doc:`installing a binary package of your choosing <install_binaries>`. Once you have done so and if you are interested in getting our Anima scripts (scripts that use Anima tools to perform complex processing series), look for the installation steps :doc:`here <install_anima_scripts>`.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: First steps

   compile_source
   install_binaries
   install_anima_scripts

Tools documentation
-------------------

More detailed documentation on the different tools available in Anima current release is provided here. It includes documentation on the following tool categories:

* :doc:`Diffusion imaging <diffusion>`
* :doc:`Registration <registration>`
* :doc:`Segmentation <segmentation>`
* :doc:`Patient to group comparison <group_comparison>`
* :doc:`MR relaxometry <relaxometry>`
* :doc:`MR denoising <denoising>`
* :doc:`Basic image processing tools <basic_processing>`


.. toctree::
   :hidden:
   :caption: Tools documentation
   
   diffusion
   registration
   segmentation
   group_comparison
   relaxometry
   denoising
   basic_processing

Scripts documentation
---------------------

And finally we provide a documentation of the main scripts available in Anima scripts that basically use Anima to make fancy things like:

* :doc:`Atlasing <atlasing>`
* :doc:`Brain extraction <brain_extraction>`
* :doc:`Diffusion imaging <diffusion_scripts>`
* :doc:`Relaxometry <relaxometry_scripts>`

.. toctree::
   :hidden:
   :caption: Scripts documentation
   
   atlasing
   brain_extraction
   diffusion_scripts
   relaxometry_scripts

Citing Anima or Anima scripts
-----------------------------

So far, there is no white paper on Anima or Anima scripts. If you are using Anima for your research, please take the time to:

* cite the relevant papers referenced in the different pages of this documentation
* add a sentence in your acknowledgments section or footnote in the methods part of your paper referencing the `RRID <https://rrid.org>`_ of Anima\: RRID\:SCR_017017 or Anima scripts\: RRID\:SCR_017072 (should be written like this to ensure it is recognized by search bots)

Do not hesitate to contact us if you have any doubts anima-users@inria.fr. Thanks !

Licensing
---------

Anima is licensed under an Aferro GPL v3 license. More details are provided in the :doc:`licensing page <license>`.

.. toctree::
   :hidden:
   :caption: Licensing
   
   license
