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
* Brain masking using the :doc:`brain extraction script <segmentation_scripts>`
* DTI estimation using the **animaDTIEstimator** tool

Some of these steps may be discarded using the ``--no-\*`` options available in the script help. The brain masking step is performed either on a provided T1 image or on the DWI first sub-volume. For distortion correction, the reversed PED image has to be provided with the ``-r`` option and the direction of the PED with the ``-d`` option. If no reversed PED image, and the T1 is, distortion will be corrected by a simple B0 to T1 non linear registration.

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

Multi-compartment along fibers analysis
---------------------------------------

We have introduced in [5] a new two parts pipeline for patient-specific, along fibers diffusion properties analysis. It seeks differences along major fiber tracts of diffusion properties between a patient MCM image and a set of controls MCM images, all registered on a common template reference for comparison. The pipeline operates in two parts:

* an offline part, done once and for all patients, that constructs from a set of controls DWIs and T1 images an MCM atlas with major fiber tracts enriched with MCM derived information
* an online part that starts from a T1 and DWI image of a patient, preprocesses them, registers them on the atlas and compares MCM enriched fibers of the patient to the ones of the controls

Fiber atlas construction
^^^^^^^^^^^^^^^^^^^^^^^^

This part of the pipeline requires that you install the TractSeg python package. For more information, please refer `here <https://github.com/MIC-DKFZ/TractSeg>`_. Once this installation is done, you will have to organize your controls data so that the actual scripts are able to use it. The data structure is the following

::

	Main_folder
	+-- T1_Data
	|  +-- PrefixT1_1.nii.gz
	|  +-- PrefixT1_2.nii.gz
	|  +-- PrefixT1_3.nii.gz
	|  .
	|  .
	|  .
	+-- Diffusion_Data
	|  +-- Diffusion_Data_1.nii.gz
	|  +-- Diffusion_Data_1.bval
	|  +-- Diffusion_Data_1.bvec
	|  +-- Diffusion_Data_1_reversed_b0.nii.gz
	|  +-- Dicom_1/ (optional)
	|  +-- Diffusion_Data_2.nii.gz
	|  +-- Diffusion_Data_2.bval
	|  +-- Diffusion_Data_2.bvec
	|  +-- Diffusion_Data_2_reversed_b0.nii.gz
	|  +-- Dicom_2/ (optional)
	|  .
	|  .
	|  .

This structure is an example with everything in. In practice, you do not need Dicom folders of the DWI data, but if you have them and have a Siemens scanner, the gradient vectors will be extracted from them to suit Anima needs. The reversed B0 has to have this specific suffis, but is also not required. It will just work better if you have them, otherwise, the T1 image will be used for distortion recovery. Once the data is structured, you can go on and run the first script (the preparation script). while being in the main folder, simply type:

.. code-block:: sh
	
	~/Anima-Scripts-Public/diffusion/mcm_fiber_atlas_comparison/animaSubjectsMCMFiberPreparation.py -n 50 -i Diffusion_Data/Diffusion_Data -d Diffusion_Data/Dicom -t T1_Data/PrefixT1

The ``-n`` option controls the number of subjects. As mentioned above, the ``-d`` option may be removed if you are confident with your gradient vectors. This script will do the following:

* preprocess each DWI using diffusion pre-processing script
* estimate MCM from smooth DWIs
* run tractseg (registration to MNI, run tractseg, get begin and end regions, merge them and put them back on patient)
* put all data in structure for atlas creation and post-processing

This script will give the following folders as output:

* Tensors: tensors for atlas creation
* Preprocessed_DWI: DWI brain masks and corrected DWI using the preprocessing scripts
* Tracts_Masks: masks for tractography from tractseg
* MCM: MCM estimations from DWI

With this preprocessed data, we can now move on to the second offline step: create a template from the control subjects. For that, create a sub-folder Atlas in your main folder and create links to the four created folders in it:

.. code-block:: sh
	
	cd Main_folder
	mkdir Atlas
	cd Atlas
	ln -s ../Tensors Tensors
	ln -s ../MCM MCM
	ln -s ../Tracts_Masks Tracts_Masks
	ln -s ../Preprocessed_DWI DWI

Run then a DTI atlas creation using the :doc:`atlasing scripts of Anima <atlasing>` in the Atlas folder with 10 iterations of atlas creation. This will build a DTI atlas from the data and output a residualDir folder to put back inside the Atlas folder as well as a averageDTI10.nrrd file. The last step of fiber atlas construction is to perform tractography on the atlas and enrich the fibers with information from the MCM models. This is done with the **animaAtlasTractsExtraction.py** script while still being in the :

.. code-block:: sh
	
	~/Anima-Scripts-Public/diffusion/mcm_fiber_atlas_comparison/animaAtlasTractsExtraction.py -i Tensors/DTI -a averageDTI10.nrrd -n 22 -m MCM/MCM_avg --mask-images-prefix DWI/DWI_BrainMask -t Tracts_Masks


This script will provide several outputs including Atlas_Tracts (raw atlas tracts on the atlas) and Augmented_Atlas_Tracts (atlas tracts with information from eachcontrol subject extracted from MCM).

Patient to atlas comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We are now ready to perform the analysis of patients using the fiber atlas constructed in the previous section. No real organization is required here, however as we may be processing several patients in a row, it is always good to get organized. Let us assume the data is organized in a patients_folder somewhere. This folder contains a DWI folder organized similarly to before and a T1 folder for anatomical data. The patient to atlas comparison may be run as follows:

.. code-block:: sh
	
	~/Anima-Scripts-Public/diffusion/mcm_fiber_atlas_comparison/animaPatientToAtlasEvaluation.py -n 50 -i DWI/Diffusion_Data_1.nii.gz -d DWI/Dicom_1 -t T1/T13D_1.nii.gz -a <Controls_Main_Folder>/Atlas/averageDTI10.nrrd -r <Controls_Main_Folder>/Atlas/Atlas_Tracts --tracts-folder <Controls_Main_Folder>/Atlas/Augmented_Atlas_Tracts

Again, the ``-d`` option is not mandatory. This will perform the following tasks using the atlas in the ``<Controls_Main_Folder>`` path:

* pre-process the patient image as in subjects preparation
* register DTI patient image on atlas DTI image, apply transformation to MCM image
* for each tract, augment raw atlas tracts with data from patient
* for each tract, compute pairwise tests along tracts between patient and controls
* Finally, compute disease burden scores for each tract

The results of all these operations will be available in two folders: **Patients_Disease_Scores** for the disease scores in csv format, **Patients_Augmented_Tracts** for the along tracts abnormalities (``*_FDR.fds``).

References
----------

1. Renaud Hédouin, Olivier Commowick, Elise Bannier, Benoit Scherrer, Maxime Taquet, Simon Warfield, Christian Barillot. *Block-Matching Distortion Correction of Echo-Planar Images With Opposite Phase Encoding Directions*. IEEE Transactions on Medical Imaging, in press available online, 2017.
2. Nicolas Wiest-Daesslé, Sylvain Prima, Pierrick Coupé, Sean Patrick Morrissey, Christian Barillot. *Rician noise removal by non-Local Means filtering for low signal-to-noise ratio MRI: applications to DT-MRI*. 11th International Conference on Medical Image Computing and Computer-Assisted Intervention. Springer, 5242 (Pt 2), pp.171-179, 2008.
3. Aymeric Stamm, Olivier Commowick, Simon K. Warfield, Simone Vantini. *Comprehensive Maximum Likelihood Estimation of Diffusion Compartment Models Towards Reliable Mapping of Brain Microstructure*. 19th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2016.
4. Aymeric Stamm, Olivier Commowick, Patrick Pérez, Christian Barillot. *Fast Identification of Optimal Fascicle Configurations from Standard Clinical Diffusion MRI Using Akaike Information Criterion*. IEEE International Symposium on Biomedical Imaging, 2014.
5. O\. Commowick, R\. Hédouin, C\. Laurent, J\.-C\. Ferré. *Patient specific tracts-based analysis of diffusion compartment models: application to multiple sclerosis patients with acute optic neuritis*. ISMRM 2021.
