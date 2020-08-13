Atlasing scripts
================

This section describes atlasing scripts of Anima scripts. These scripts are special in Anima scripts as they use shell and the `OAR <http://oar.imag.fr>`_ job submission system.

Overall goal
------------

Two scripts are available here that have to be run on a OAR cluster to work properly. They follow a modified Guimond et al. [1] method to compute an average unbiased atlas from respectively anatomical or diffusion tensor images.

The method proceeds iteratively by registering the images onto a current reference to compute the reference at the next iteration. The main modifications with respect to the original Guimond et al. method is to use the log-Euclidean framework on diffeomorphisms [2] to compute the average transformations. 

The two scripts then differ on the registration method used: anatomical image non linear registration [3] or DTI non linear registration [4].

Data structure for atlas construction
-------------------------------------

The datasets are supposed to be named correctly and placed in the right folders so that the atlasing scripts will find them. Data should be organized like this:

::

	Main folder
	+-- Images
	|  +-- Prefix_1.nii.gz
	|  +-- Prefix_2.nii.gz
	|  +-- Prefix_3.nii.gz
	|  .
	|  .
	|  .
	+-- Masks
	|  +-- Mask_1.nii.gz
	|  +-- Mask_2.nii.gz
	|  +-- Mask_3.nii.gz
	|  .
	|  .
	|  .

The main folder can be located anywhere, the prefix of the images names (i.e. **Prefix**) and their folder name (i.e. **Images**) are up to the user (but should not include spaces, special characters or accentuated characters by precaution). The images folder contains the images from which the atlas will be built (anatomical or DTI depending on the script), all suffixed from 1 to N without losing contiguity. 

The Masks folder and the subsequent names should be as shown here. They concern however only the anatomical images script for which they provide the mask of interest for each image that will be considered for average image construction.

Atlasing scripts
----------------

With this folder structure, two scripts are available that follow the same idea: **animaBuildAnatomicalAtlas.py** for anatomical atlas creation and animaBuildDTIAtlas.py** for DTI atlas creation. Both should be run from an OAR frontend machine as they start many jobs (N per iteration where N is the number of images, times a total of M iterations). 

To run them, use one of the commands that follow after going to the main folder in the terminal (depending if you are making a DTI atlas or anatomical one):

.. prompt:: bash $

	~/Anima-Scripts-Public/atlasing/anatomical/animaBuildAnatomicalAtlas.py -p Images/Prefix -n <N> -i <M> -c <nCores>
	~/Anima-Scripts-Public/atlasing/dti/animaBuildDTIAtlas.py -p Images/Prefix -n <N> -i <M> -c <nCores>

This runs an atlas construction over N images, performing M iterations of atlas creation and asking for nCores cores on each computer of the cluster for each job. 

Results are provided in the main folder in the form of an average image (**averageForm<M>.nii.gz** or **averageDTI<M>.nii.gz** depending on the script). Transformations bringing each individual image on the atlas can be found in the **residualDir** subfolder added by the script:

* `Prefix_<i>_linear_tr.txt`: affine transform from image i toward **averageForm<M-1>.nii.gz**
* `Prefix_<i>_nonlinear_tr.nrrd`: non linear transform from affinely transformed image i toward **averageForm<M-1>.nii.gz**
* `sumNonlinear_inv_tr.nrrd`: average transformation whose inverse should be applied to get transformed images onto **averageForm<M>.nii.gz**

Online atlasing (iterative centroid)
------------------------------------

Online atlasing as explained in [6] works in a different way as the previous scripts. Assuming you hav an atlas output products from the previous scripts, the iterative centroid algorithm adds additional images to form an atlas without having to create it from the ground up. Its call works in the same way as the previous ones:

.. prompt:: bash $

	~/Anima-Scripts-Public/atlasing/anatomical_iterative_centroid/animaBuildAnatomicalICAtlas.py -p Images/Prefix -s <startPoint> -n <N> -c <nCores>

The `-s` option gives the number of images already in the atlas. The script will then use the original atlas and additional images up to N to compute the new atlas. It assumes the same data structure as the previous scripts

Longitudinal atlasing
---------------------

A longitudinal atlas is constituted of a set of 3D atlases each representing an average model at a given age. In this case, the 3D atlases are computed up to a rigid transformation accounting both for global and local transformation in the unbiasing step.

We provide in the following sub-sections the necessary scripts to compute longitudinal atlases as detailed in [5], strongly based on the previous scripts. A preliminary weighting step is necessary such that, in the creation of an atlas at time *t*, more importance is given to subjects closer in age to *t*.

Weighting
^^^^^^^^^

This script takes as an input all of the following:

* List of all images paths in a txt file (i.e. image.txt). One path per line
* List of associated ages in an other txt file (i.e. age.txt)
* List of desired atlas ages in an other txt file (i.e. atlasAge.txt)

.. prompt:: bash $

	~/Anima-Scripts-Public/atlasing/longitudinal_preparation/animaComputeLongitudinalAtlasWeights.py -a age.txt -i image.txt -o outputFolder -n n -A atlasAge.txt -p Images/Prefix

This script call runs the preparation of everything needed to compute a 4D atlas composed of a set of 3D atlases representatives of ages contained in atlasAge.txt, made using about *n* subjects. It also prepares the folder structure to compute the different atlases in the following step.

::

	output folder
	+-- atlas_1
	.
	.
	.
	+-- atlas_i
	|  +-- Images
	|  |  +-- Prefix_1.nii.gz
	|  |  .
	|  |  .
	|  |  .
	|  +-- weights.txt
	.
	.
	.

Atlasing scripts with longitudinal parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the preparation step, to compute each sub-atlas i, simply run one of the atlasing scripts with the appropriate options:

.. prompt:: bash $

	cd outputFolder/atlas_i
	~/Anima-Scripts-Public/atlasing/anatomical/animaBuildAnatomicalAtlas.py -p Images/Prefix -n <N> -i <M> -c <nCores> --rigid -w weights.txt -b 2

In the previous line, **animaBuildAnatomicalAtlas** may be replaced by **animaBuildDTIAtlas** to compute a DTI longitudinal atlas.

References
----------

1. A. Guimond, J. Meunier, J.P. Thirion. *Average brain models: A convergence study*, Computer Vision and Image Understanding, 77(2):192-210, 2000.
2. V. Arsigny, O. Commowick, X. Pennec, N. Ayache. *A Log-Euclidean Framework for Statistics on Diffeomorphisms*, 9th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 924-931, 2006.
3. Olivier Commowick, Nicolas Wiest-Daesslé, Sylvain Prima. *Automated diffeomorphic registration of anatomical structures with rigid parts: application to dynamic cervical MRI*. 15th International Conference on Medical Image Computing and Computer Assisted Intervention, pp.163-70, 2012.
4. Ralph Suarez, Olivier Commowick, Sanjay Prabhu, Simon K. Warfield. *Automated delineation of white matter fiber tracts with a multiple region-of-interest approach*. NeuroImage, 59 (4), pp.3690-3700, 2012.
5. Antoine Legouhy, Olivier Commowick, François Rousseau, Christian Barillot. *Unbiased Longitudinal Brain Atlas Creation Using Robust Linear Registration and Log-Euclidean Framework for Diffeomorphisms*, International Symposium on Biomedical Imaging, 2019.
6. Antoine Legouhy, Olivier Commowick, François Rousseau, Christian Barillot.  *Online Atlasing Using an Iterative Centroid*, MICCAI, 2019.