Segmentation scripts
====================

This section describes Anima scripts for segmentation. For now, we have atlas-based intracranial extraction and multi-atlas segmentation.

Atlas-based intracranial mask extraction
----------------------------------------

This requires an additional atlas image for which the intracranial mask is known. These data are provided in the Anima scripts data.

The script simply takes as an input the image or images to be brain-extracted and performs a series of registrations to bring the atlas on the image. Finally the brain mask is applied to the image. For each brain-extracted image, two images are produced: ``image_brainMask.nrrd`` and ``image_masked.nrrd``. The first one is the binary mask of the brain, and the second one the image with only the brain kept. T1 images are the ones that should work best as the atlas proposed is in the T1 modality but since the registration algorithm uses an adapted similarity measure, other MRI modalities should be fine as well as long as their field of view is similar to the one of the atlas.

*Example:*

.. code-block:: sh
	
	~/Anima-Scripts-Public/brain_extraction/animaAtlasBasedBrainExtraction.py -i T1Image.nrrd

You may also now specify a folder with the ``-a`` option that contains another atlas than the one used by default in the script. The folder must contain three files to work:

* ``Reference_T1.nrrd``: the atlas T1w image, not brain masked
* ``Reference_T1_masked.nrrd``: the brain masked atlas T1w image
* ``BrainMask.nrrd``: the atlas brain mask

Multi-atlas segmentation
------------------------

As for the :doc:`atlasing <atlasing>` scripts, this script requires to be run on a cluster with an `OAR <http://oar.imag.fr>`_ scheduler. It provides a basic multi-atlas segmentation, i.e. does the following:

* registers a set of images with known segmentations (atlases) onto an image to be delineated
* applies transformations to known segmentations
* fuses the transformed segmentations using majority voting (**animaMajorityLabelVoting** in :doc:`segmentation tools <segmentation>`)

Several options are available:

* ``-i``: anatomical image to be delineated
* ``-a``: list of atlas anatomical images i.e. a text file with a file name on each line
* ``-s``: list of corresponding atlas segmentations i.e. a text file with a segmentation name on each line
* ``-o``: output label image for the input anatomical image
* ``-c``: optional number of cores for each job on the cluster

*Example:*

.. code-block:: sh
	
	~/Anima-Scripts-Public/multi_atlas_segmentation/animaMultiAtlasSegmentation.py -i T13D.nrrd -a listAtlasImages.txt -s listAtlasSegmentations.txt -o T13D_segmented.nrrd

Tissues classification
----------------------

This script performs the task of tissues classification of a brain image (possibly with several modalities) using an external atlas of tissue probabilities. It basically performs the regsitrations needed to bring the atlas onto the brain image to classify and then run tissues classification using :doc:`animaTissuesEMSegmentation <segmentation>`. The script name is **animaAtlasEMTissuesSegmentation** and can be used as follows.

*Example:*

.. code-block:: sh
	
	~/Anima-Scripts-Public/em_segmentation/animaAtlasEMTissuesSegmentation.py -i T1Image.nrrd -i T2Image.nrrd -m brain_mask.nrrd -o output_classification.nrrd

Among the options used here, are the following:

* ``-i``: anatomical image(s) of the brain to be segmented. The images need not be registered to each other
* ``-m``: brain mask of the first input
* ``-o``: output file name for the brain classification
