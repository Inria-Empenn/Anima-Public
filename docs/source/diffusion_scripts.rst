Diffusion imaging scripts
=========================

This section describes diffusion imaging preprocessing and model estimation scripts of Anima scripts. 

Diffusion images preprocessing script
-------------------------------------

This script combines several preprocessing steps to prepare data with the following steps (in that order):

* if -D option is given, correct gradients to match requirements from Anima (in real space i.e. the scanner space)
* Eddy current correction and motion correction using the experimental tool from Anima **animaEddyCurrentCorrection**
* distortion correction using the method proposed in [1]
* reorientation of the DWI volume to be axial on the z-axis and have no reversed axis
* Denoising using the NL-Means method [2]
* Brain masking using the :doc:`brain extraction script <brain_extraction>`
* DTI estimation using the **animaDTIEstimator** tool

Some of these steps may be discarded using the --no-\* options available in the script help. The brain masking step is performed either on a provided T1 image or on the DWI first sub-volume. For distortion correction, the reversed PED image has to be provided with the -r option and the direction of the PED with the -d option. If no reversed PED image, and the T1 is, distortion will be corrected by a simple B0 to T1 non linear registration.

Warning: gradients reworking is known to work only with Siemens acquisitions, not tested on other scanners.

*Example:*

.. code-block:: sh
	
	~/Anima-Scripts-Public/diffusion/animaDiffusionImagePreprocessing.py -b Diff.bval -D Dicom/* -r B0_PA.nii.gz -d 1 -t T1.nii.gz -i Diff.nii.gz

Results of this scripts are:

* **Diff_preprocessed.nrrd**: preprocessed diffusion 4D image
* **Diff_preprocessed.bvec**: preprocessed diffusion gradient vectors
* **Diff_brainMask.nrrd**: diffusion brain mask
* **Diff_Tensors.nrrd**: estimated tensors
* **Diff_Tensors_B0.nrrd**: estimated tensors B0 image
* **Diff_Tensors_NoiseVariance.nrrd**: estimated noise variance from the tensor estimation (actually includes noise and model underfit error)

Multi-compartment models estimation script
------------------------------------------

This script is more experimental as the multi-compartment model estimation tool in Anima [3] is currently undergoing quite some changes. However, it works well and we have a script that uses preferably the results of the diffusion preprocessing script as an input.

It mainly has two modes: one for HCP-like datasets that have high quality acquisitions and enable less constrained models, and one for regular clinical data with less estimation demanding models. Several options are still available:

* **-t**: choose the model type (as explained in :doc:`diffusion documentation <diffusion>`
* **-n**: maximal number of anisotropic compartments (ideally choose a number in between 1 and 3)
* **--hcp**: add some additional compartments like stationary water to handle specific aspects of `HCP data <https://www.humanconnectome.org>`_
* **--no-model-simplification** and **-S** control model selection / averaging option. If none are set, the script will by default follow the method proposed in [4]: compute all models from 0 to N compartments and "average" them according to their likelihood (according to the AIC criterion). If **-S** is set, a faster model selection done at the stage of the stick model estimation is performed. If **--no-model-simplification** is set, a pure N model estimation is done without model selection.

*Example:*

.. code-block:: sh

	~/Anima-Scripts-Public/diffusion/animaMultiCompartmentModelEstimation.py -t tensor -n 3 -i Diff_preprocessed.nrrd -g Diff_preprocessed.bvec -b Diff.bval -m Diff_brainMask.nrrd

References
----------

1. Renaud Hédouin, Olivier Commowick, Elise Bannier, Benoit Scherrer, Maxime Taquet, Simon Warfield, Christian Barillot. *Block-Matching Distortion Correction of Echo-Planar Images With Opposite Phase Encoding Directions*. IEEE Transactions on Medical Imaging, in press available online, 2017.
2. Nicolas Wiest-Daesslé, Sylvain Prima, Pierrick Coupé, Sean Patrick Morrissey, Christian Barillot. *Rician noise removal by non-Local Means filtering for low signal-to-noise ratio MRI: applications to DT-MRI*. 11th International Conference on Medical Image Computing and Computer-Assisted Intervention. Springer, 5242 (Pt 2), pp.171-179, 2008.
3. Aymeric Stamm, Olivier Commowick, Simon K. Warfield, Simone Vantini. *Comprehensive Maximum Likelihood Estimation of Diffusion Compartment Models Towards Reliable Mapping of Brain Microstructure*. 19th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2016.
4. Aymeric Stamm, Olivier Commowick, Patrick Pérez, Christian Barillot. *Fast Identification of Optimal Fascicle Configurations from Standard Clinical Diffusion MRI Using Akaike Information Criterion*. IEEE International Symposium on Biomedical Imaging, 2014.
