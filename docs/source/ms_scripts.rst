Multiple sclerosis scripts
==========================

This section describes Anima scripts created to segment lesions from brain MRI images in multiple sclerosis (MS) patients.

MSSEG 2016 preprocessing script
-------------------------------

This script (**animaMSExamPreparationMSSEG2016.py**) combines several processing steps to preprocess a single time point containing T1, T2/PD, FLAIR and Gd T1 images using the same script that was used in the MSSEG 2016 challenge [1]. It relies on a precomputed brain mask that may be obtained from `volBrain <https://www.volbrain.upv.es>`_ (as for the challenge) or :doc:`the brain extraction scripts of Anima <segmentation_scripts>`. The script performs denoising, bias correction, brain masking and registration to the FLAIR image as explained in [1].

MSSEG-2 longitudinal preprocessing script
-----------------------------------------

This second script comes from the `second challenge <https://portal.fli-iam.irisa.fr/msseg-2/>`_ on new MS lesions segmentation that took place at MICCAI 2021. It is named **animaMSLongitudinalPreprocessing.py** and takes as an input an organized folder with FLAIR images at two time points (already rigidly registered) and, if available, of some ground truths. It then performs basic preprocessing using Anima tools, on the two images. This script was provided to the challengers at the MSSEG-2 challenge to help them participate. The following steps are performed:

* brain mask extraction on the two images and fusion
* brain masking of the images
* N4 bias correction
* if a template FLAIR image is given, Nyul normalization of the images to the template

References
----------

1. Olivier Commowick, Michaël Kain, Romain Casey, Roxana Ameli, Jean-Christophe Ferré, Anne Kerbrat, Thomas Tourdias, Frédéric Cervenansky, Sorina Camarasu-Pop, Tristan Glatard, Sandra Vukusic, Gilles Edan, Christian Barillot, Michel Dojat, Francois Cotton. *Multiple sclerosis lesions segmentation from multiple experts: The MICCAI 2016 challenge dataset*. NeuroImage, 244, pp.1-8, 2021.
