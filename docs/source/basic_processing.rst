Basic image processing tools
============================

Image arithmetic
----------------

**animaImageArithmetic** allows you to add, subtract, multiply or divide images either by a constant or another image. It relies on ITK filters.

*Example:* this computes Image1 - Image2 into the result image Diff.nii.gz

.. code-block:: sh

	animaImageArithmetic -i Image1.nii.gz -s Image2.nii.gz -o Diff.nii.gz

ROI intensities statistics
--------------------------

**animaROIIntensitiesStats** provides statistics (mean, max, median, variance) in regions of interests on an image or a list of images in a text file.

*Example:* this computes the median in each label of roi_image.nii.gz of the input image image.nii.gz. The output is given in a csv file data.csv

.. code-block:: sh

	animaROIIntensitiesStats -i image.nii.gz -r roi_image.nii.gz -o data.csv -s median

Image averaging
---------------

**animaAverageImages** provides a tool to average multiple images (vector or scalar) into a single one. 

*Example:* this computes the average from all images listed in listFiles.txt. listMasks.txt specifies a list of masks of interest applicable for each image.

.. code-block:: sh
 
	animaAverageImages -i listFiles.txt -m listMasks.txt -o outputImage.nrrd 

Image collapse
--------------

**animaCollapseImage** is a simple tool for converting a vector image (e.g. 3D vector file) into a (N+1)D image (e.g. a 4D scalar image).

Images concatenation
--------------------

**animaConcatenateImages** takes several images and concatenates them into an (N+1)D image. Options allow to set the origin and spacing of the new dimension, as well as an optional base (N+1)D volume to which the additional images will be concatenated.

*Example:* this will concatenate Base, Addon1 and Addon2 into a single (N+1)D image on the last dimension.

.. code-block:: sh

	animaConcatenateImages -i Addon1.nii.gz -i Addon2.nii.gz -b Base.nii.gz -o Data.nii.gz 

Image resolution changer
------------------------

**animaImageResolutionChanger** provides a way to change the image resolution (in the sense voxel resolution in mm). It automatically recomputes the origin and adapted image size from the old and new resolutions. The tool works for 3D (scalar or vector image), and 4D scalar images. For 4D images, the resolution change is made on  sub 3D images.

*Example:* this resamples the input image so that its voxel resolution is 2x2x2 mm.

.. code-block:: sh

	animaImageResolutionChanger -i Image.nrrd -o Image_Result.nrrd -x 2 -y 2 -z 2

Image vectorization
-------------------

**animaVectorizeImages** takes several ND scalar volumes and outputs a single ND vector image. It can also take a text file listing several volumes as the input.

*Example:* this will create a vector image with vector dimension 2, where the first dimension will contain Image1 and the second Image2.

.. code-block:: sh

	animaVectorizeImages -i Image1.nii.gz -i Image2.nii.gz -o VectImage.nii.gz

Connected components
--------------------

**animaConnectedComponents** computes connected components from a binary image. Optionally, it may take a minimal component size, under which the components will be removed.

Image conversion
----------------

**animaConvertImage** serves several roles: it can display information on an image, reorient it and write it to another file format.

*Example:* this reorients the image in the coronal plane, displays its information (``-I``) and saves the result in a NRRD file.

.. code-block:: sh

	animaConvertImage -i Image.nii.gz -I -R CORONAL -o Image_Coronal.nrrd

Shapes format conversion
------------------------

**animaConvertShapes** allows you to convert shapes (fibers, surfaces, etc.) between file formats supported by Anima (vtk, vtp, fds, csv).

*Example:* this converts a VTK ascii file to a VTP file.

.. code-block:: sh

	animaConvertShapes -i Shape.vtk -o Shape.vtp

Image cropping
--------------

**animaCropImage** extracts a sub-volume from an image using the ITK ExtractImageFilter. The lower case arguments(x<xindex>, y<yindex>, z<zindex>, t<tindex>) are the starting indexes of the input region to keep. The default value is 0. The upper case arguments (X<xsize>, Y<ysize>, Z<zsize>, T<tsize>) are the sizes of the input region to keep. The default value is the largest possible size given the corresponding indexes.

If you give arguments size of zero the corresponding dimension will be collapsed.

*Example:* for a 4D image 4x4x4x4, the arguments ``--xindex 1 --zindex 1 --zsize 2 --tindex 3 --tsize 0`` will result in an image 3x4x2 where the x dim corresponds to [1,2,3] of the input, y[0,3], zindex[1,2] and tindex is collapsed, only the last sequence has been kept.

Image smoothing
---------------

**animaImageSmoother** simply applies Gaussian smoothing with a specific sigma value to an image using Young - Van Vliet's recursive smoothing filter implemented in ITK [1].

Morphological operations
------------------------

**animaMorphologicalOperations** computes usual morphological operations (erosion, dilation, opening, closure), with a specified radius expressed in millimeters. 

References
----------

1. Irina Vidal-Migallón, Olivier Commowick, Xavier Pennec, Julien Dauguet, Tom Vercauteren. *GPU & CPU implementation of Young - Van Vliet's Recursive Gaussian Smoothing Filter*. Insight Journal (ITK), 2013, pp.16
