Brain extraction
================

This section describes the brain extraction scripts of Anima scripts. 

Atlas-based intracranial mask extraction
----------------------------------------

This requires an additional atlas image for which the intracranial mask is known. These data are provided in the Anima scripts data.

The script simply takes as an input the image or images to be brain-extracted and performs a series of registrations to bring the atlas on the image. Finally the brain mask is applied to the image. For each brain-extracted image, two images are produced: **image_brainMask.nrrd** and **image_masked.nrrd**. The first one is the binary mask of the brain, and the second one the image with only the brain kept. T1 images are the ones that should work best as the atlas proposed is in the T1 modality but since the registration algorithm uses an adapted similarity measure, other MRI modalities should be fine as well as long as their field of view is similar to the one of the atlas.

*Example:*

.. code-block:: sh
	
	~/Anima-Scripts-Public/brain_extraction/animaAtlasBasedBrainExtraction.py T1Image.nrrd
