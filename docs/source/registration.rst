Registration tools
==================

This part describes registration features in Anima. It includes linear and non linear registration (of anatomical and DTI images), as well as EPI distortion correction and tools for applying series of transformations to images. One note on transformations that come out of these tools: they are all in real coordinates (as opposed to voxel coordinates).

Linear registration
-------------------

Anima provides linear registration (**animaPyramidalBMRegistration**) using a block-matching algorithm as presented in [1,2]. It may use BOBYQA or exhaustive optimization for finding matches (``--opt`` parameter). Transformations between blocks may be translations, rigid or affine (experimental) transformations (``-t`` parameter). The result global transform (writable in any ITK supported format) may be a translation, rigid or affine transformation. Both this executable and the next one use a local squared correlation coefficient as the default similarity measure between blocks (``--metric`` parameter). Registration may be initialized from a previously computed transformation (``-i`` option) or from intensity based PCA translation or rigid transformation (``-I`` option). In addition, we have introduced new constrained affine transformations [13] (``--ot`` option) enabling the quantification of growth factors of e.g. the brain on specific chosen directions, and the extraction of the nearest rigid or similarity transformation from an affine registration (``--out-rig`` and ``--out-sim`` options) [14].

*Example:* this registers Floating.nii.gz on Reference.nii.gz using a global affine transform (``--ot`` parameter), on a multi-resolution pyramid of 4 levels (``-p`` parameter), but stopping at the one before last level (``-l`` parameter). For matching, it uses only those blocks where the standard deviation is above 10.

.. code-block:: sh

	animaPyramidalBMRegistration -r Reference.nii.gz -m Floating.nii.gz -o Floating_On_Ref.nii.gz -O transform_aff.txt -s 10 -p 4 -l 1 --ot 2

Non linear registration
-----------------------

Anatomical registration
^^^^^^^^^^^^^^^^^^^^^^^

Non linear registration of anatomical images (**animaDenseSVFBMRegistration**) shares many parameters with linear registration. It however estimates a stationary velocity field between two images [3], modeling a dense non linear  transformation between the images. The output is therefore a vector image containing the SVF. This executable implements the method proposed in [4]. It requires that the two input images have the same sizes and orientations. It is usually highly recommended to perform non linear registration after global rigid/affine registration.

*Example:* this registers Floating_aff.nii.gz on Reference.nii.gz on a multi-resolution pyramid of 3 levels (``-p`` parameter), stopping at full resolution (``-l`` parameter). 

.. code-block:: sh

	animaDenseSVFBMRegistration -r Reference.nii.gz -m Floating_aff.nii.gz -o Floating_On_Ref.nii.gz -O transform_nl.nii.gz -p 3 -l 0 

DTI registration
^^^^^^^^^^^^^^^^

Non linear DTI registration (**animaDenseTensorSVFBMRegistration**) implements the same algorithms as in scalar images registration with finite strain tensor reorientation or preservation of principal direction (PPD). Parameters are the same except the similarity metrics which implement those proposed in [4] (tensor oriented generalized correlation). Note that the block variance may be much smaller for tensors and it is advised to set the minimal variance (``-s``) to 0 by default.

MCM registration
^^^^^^^^^^^^^^^^

Non linear MCM (multi-compartment models) registration (**animaDenseMCMSVFBMRegistration**) implements the same algorithms as in scalar images registration with finite strain tensor reorientation or preservation of principal direction (PPD). Parameters are the same except similarity metrics which implement those proposed in [12] (MCM SSD and correlation surrogate). As for tensors the block variance may be much smaller than for scalar images and it is advised to set the minimal variance (``-s``) to 0 by default.

EPI artifacts correction
------------------------

Eddy current distortion correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Anima includes a simple, yet experimental, tool for Eddy current distortion correction **animaEddyCurrentCorrection**. It performs by registering linearly then non linearly each sub-volume of the EPI series to the first sub-volume. Each non linear transformation is computed only in the phase encoding direction as suggested by other publications. Please note that the Eddy current distortion correction tool does not include susceptibility distortion correction that is accounted for by the distortion correction tool in the next section.

*Example:* this will perform Eddy current distortion correction on DiffVolume4D.nrrd in the direction specified by ``-d`` (``-d 1`` denotes the Y direction in voxel coordinates).

.. code-block:: sh

	animaEddyCurrentCorrection -i DiffVolume4D.nrrd -o DiffVolume_Corrected.nrrd -d 1 

Susceptibility distortion correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Anima proposes two tools for EPI distortion correction of images: 

* **animaDistortionCorrection** which implements the method proposed by Voss et al. [5]
* **animaBMDistortionCorrection** which implements the symmetric block-matching method proposed in [6]. It uses the same kind of parameters as previous block-matching registration algorithms.

From our experience, we found out it is best to use both for correcting for distortion: **animaDistortionCorrection** being used as the initialization for **animaBMDistortionCorrection**. Both tools use two 3D images with reversed phase encoding directions as their input, the direction that was used in the acquisition may be specified using the ``-d`` parameter (0 = X axis, 1 = Y axis, 2 = Z axis). The transformation found may then be applied to a whole 4D volume using our resampling tools.

*Example:* this first computes an initial correction field using Voss et al. method, with a smoothing sigma of 2 voxels (``-s`` parameter). Then it starts from the initial transformation Init_Correction.nii.gz and computes a more precise BM_Correction.nii.gz dense transformation (note that contrary to previous algorithms, these transformations are not SVFs but dense displacement fields). It also outputs the average of AP and PA images after correction into BM_Correction.nii.gz.

.. code-block:: sh

	animaDistortionCorrection -f AP_Image.nii.gz -b PA_Image.nii.gz -o Init_Correction.nii.gz -s 2
	animaBMDistortionCorrection -f AP_Image.nii.gz -b PA_Image.nii.gz -i Init_Correction.nii.gz -o BM_Corrected_Image.nii.gz -O BM_Correction.nii.gz

Symmetry plane computation and constrained registration
-------------------------------------------------------

In addition to traditional registration, we provide tools to compute and use the inter-hemispheric symmetry plane of an image [7,8]. This is based on two tools:

* **animaSymmetryPlane** [7] computes the symmetry transformation of an image (about its inter-hemispheric plane) and outputs both that transform (``-O`` parameter) and a transformation that brings the image on its symmetry plane (``--out-realign-trsf``)
* **animaSymmetryConstrainedRegistration** implements constrained global rigid registration [8] utilizing two input symmetry plane transforms to restrict the search space.

*Example:* If one wants to register two images A.nii.gz and B.nii.gz, three steps will be necessary: realign A on its symmetry plane, realign B on its symmetry plane, and use both transformations as inputs to make a constrained registration of A and B. The output transformation brings the original B on the original A with a rigid transformation. The ``-F`` option activates a faster constrained registration but which may lose a little accuracy (see [8]).

.. code-block:: sh

	animaSymmetryPlane -i A.nii.gz -o A_realign.nii.gz --out-realign-trsf A_sym.txt
	animaSymmetryPlane -i B.nii.gz -o B_realign.nii.gz --out-realign-trsf B_sym.txt
	animaSymmetryConstrainedRegistration -r A.nii.gz -m B.nii.gz --ref-sym A_sym.txt --moving-sym B_sym.txt -F -o B_on_A.nii.gz -O B_on_A_rig.txt

Transformation tools (applying, arithmetic, jacobian)
-----------------------------------------------------

EPI distortion correction
^^^^^^^^^^^^^^^^^^^^^^^^^

EPI distortion correction works in a slightly different way as other resampling tools. The tool provided is called **animaApplyDistortionCorrection**. It takes as inputs a 4D image with regular phase encoding direction (``-f`` parameter) and optionally a 4D image with reversed phase encoding direction (for better correction, ``-b`` parameter). Then, using transformations coming from the previous tools, it corrects for distortion (if ``-b`` is provided the output will be the average of the two corrected images).

*Example:* this applies the previously obtained transormation to the whole DWI volume to correct its distortion.

.. code-block:: sh

	animaApplyDistortionCorrection -f DWI_AP.nii.gz -t BM_Correction.nii.gz -o DWI_Corrected.nii.gz

Constructing series of transformations descriptions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All other transform application tools require the input transformations to be given as an XML file which describes a series of transformations. It can take several option but the simple example is the following:

.. code-block:: sh

	animaTransformSerieXmlGenerator -i transform_aff.txt -i transform_nl.nii.gz -o transforms.xml

It creates the description of the two transformations (the specified order is the order in which they will be applied).

Applying a transformation to images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Applying a transformation requires the previous XML description file. Three tools are available:

* one for scalar images - **animaApplyTransformSerie**
* one for tensor images - **animaTensorApplyTransformSerie**
* one for multi-compartment model images [10] **animaMCMApplyTransformSerie**

All tools require a geometry image to tell in which space the resampling will take place (``-g`` parameter). If the transformation series is globally linear, it may be applied to a gradient file of diffusion images. **animaApplyTransformSerie** now supports 3D and 4D images (in the latter case, the transformation is applied independently to each of the 3D sub-volumes). Diffusion model resamplers have an option to either apply finite strain re-orientation of the models or preservation of principal direction (PPD) re-orientation: ``-P`` activates PPD re-orientation, the default is finite strain.

*Example:* this applies the transforms in transforms.xml to resample Floating on Reference.

.. code-block:: sh

	animaApplyTransformSerie -i Floating.nii.gz -g Reference.nii.gz -t transforms.xml -o F_resampled.nii.gz

Applying a transformation to fibers or meshes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Anima also comes with a tool to apply transformations obtained from image registration to meshes or fibers. It accepts vtk, vtp and fds (our fiber format for `medInria <https://med.inria.fr>`_) files. This tool, named **animaFibersApplyTransformSerie**, works in the same way as animaApplyTransformSerie. The two main differences are the following:

* the input transformation is inverted by default as image transformations are encoded in Anima as the inverse of the underlying space transformation. This way, animaApplyTransformSerie and animaFibersApplyTransformSerie are similar in their uses. Use the ``-I`` option to invert the transformation series if necessary.
* There is no need for a geometry as this is specific to images

Computing the Jacobian of a transformation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A tool to compute the Jacobian or its determinant of a displacement field transformation is provided with the tool **animaDisplacementFieldJacobian**. If may handle SVF transformations using the ``-S`` option. More options for this tool are provided when using the ``-h`` option.

*Example:* this computes the Jacobian matrix of the input SVF after its exponentiation. The Jacobian matrix is stored as a 9 component vector image stored in row first. 

.. code-block:: sh

	animaDisplacementFieldJacobian -i inputField.nrrd -S -o dispFieldJacDeterminant.nrrd

Linear transformations arithmetic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide a tool named **animaLinearTransformArithmetic** to compose and perform log-Euclidean operations on linear transformations as proposed by Arisgny et al. [9]. The tool proposes regular composition (``-c``), addition (``-a``), subtraction (``-s``), multiplication by a constant (``-M``), division by a constant (-D) in the log-Euclidean framework. 

*Example:* this performs the log-Euclidean addition of the two linear input transformations (in the ITK format) in the log-Euclidean framework.

.. code-block:: sh

	animaLinearTransformArithmetic -i transform.txt -a addedTransform.txt -o outputTransform.txt

Dense field transformations arithmetic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide a tool named **animaDenseTransformArithmetic** to compose and perform log-Euclidean operations on dense field (diffeomorphic) transformations as proposed by Arisgny et al. [11]. The tool proposes regular composition or BCH approximation to the composition of SVFs in the log-Euclidean framework (``-c``). It can also take the exponential of an SVF, dense diffeomorphism logarithm is the only operation not implemented yet.

References
----------

1. Olivier Commowick, Nicolas Wiest-Daesslé, Sylvain Prima. *Block-Matching Strategies for Rigid Registration of Multimodal Medical Images*. 9th IEEE International Symposium on Biomedical Imaging (ISBI), pp. 700--703, 2012.
2. S\. Ourselin, A\. Roche, S\. Prima and N\. Ayache. *Block Matching: A General Framework to Improve Robustness of Rigid Registration of Medical Images*. Third International Conference on Medical Robotics, Imaging And Computer Assisted Surgery (MICCAI), volume 1935 of LNCS, pp. 557--566, 2000.
3. Olivier Commowick, Nicolas Wiest-Daesslé, Sylvain Prima. *Automated diffeomorphic registration of anatomical structures with rigid parts: application to dynamic cervical MRI*. 15th International Conference on Medical Image Computing and Computer Assisted Intervention, pp.163-70, 2012.
4. Ralph Suarez, Olivier Commowick, Sanjay Prabhu, Simon K. Warfield. *Automated delineation of white matter fiber tracts with a multiple region-of-interest approach*. NeuroImage, 59 (4), pp.3690-3700, 2012.
5. H.U. Voss, R. Watts, A.M. Ulugc, D. Ballona. *Fiber tracking in the cervical spine and inferior brain regions with reversed gradient diffusion tensor imaging*. Magnetic Resonance in Medicine, 24(3):231–239, 2006.
6. Renaud Hédouin, Olivier Commowick, Elise Bannier, Benoit Scherrer, Maxime Taquet, Simon Warfield, Christian Barillot. *Block-Matching Distortion Correction of Echo-Planar Images With Opposite Phase Encoding Directions*. IEEE Transactions on Medical Imaging, in press available online, 2017.
7. S\. Prima, S\. Ourselin, N\. Ayache. *Computation of the Mid-Sagittal Plane in 3D Brain Images*. IEEE Transactions on Medical Imaging, 21(2):122-138, February 2002\.
8. Sylvain Prima, Olivier Commowick. *Multimodal rigid-body registration of 3D brain images using bilateral symmetry*. Medical Imaging: Image Processing, SPIE, 8669, pp.866911, 2013.
9. V\. Arsigny, O\. Commowick, N\. Ayache, X\. Pennec. *A Fast and Log-Euclidean Polyaffine Framework for Locally Linear Registration*. Journal of Mathematical Imaging and Vision, 33(2):222-238, February 2009.
10. Renaud Hédouin, Olivier Commowick, Aymeric Stamm, Christian Barillot. *Interpolation and Averaging of Multi-Compartment Model Images*, 18th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 354-362, 2015.
11. V\. Arsigny, O\. Commowick, X\. Pennec, N\. Ayache. *A Log-Euclidean Framework for Statistics on Diffeomorphisms*, 9th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 924-931, 2006.
12. O\. Commowick, R\. Hédouin, E\. Caruyer, C\. Barillot. *L2 Similarity Metrics for Diffusion Multi-Compartment Model Images Registration*, 20th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 257-265, 2017.
13. A\. Legouhy, O\. Commowick, F\. Rousseau, C\. Barillot. *Anisotropic similarity, a constrained affine transformation: application to brain development analysis*, ISMRM, 2018.
14. A\. Legouhy, O\. Commowick, F\. Rousseau, C\. Barillot. *Unbiased Longitudinal Brain Atlas Creation Using Robust Linear Registration and Log-Euclidean Framework for Diffeomorphisms*, International Symposium on Biomedical Imaging, 2019.
