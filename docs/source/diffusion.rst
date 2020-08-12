Diffusion imaging tools
=======================

This section describes only the core estimation and tractography features of Anima. For registration of diffusion models and EPI distortions correction, please refer to the :doc:`registration page <registration>`.

Diffusion imaging
-----------------

Multiple compartment models
^^^^^^^^^^^^^^^^^^^^^^^^^^^

From version 2.0, Anima provides tools for the estimation and processing of multi-compartment models. So far, we have implemented results of our research from papers on estimation and model averaging of MCM. More works on tractography are to come in future releases.

MCM estimation
""""""""""""""

The **animaMCMEstimator** tool provides an implementation of MCM estimation from DWI data, following the maximum likelihood method described in [8]. It provides in its public version tools to estimate models including isotropic compartments (free water, isotropic restricted water, dot compartment, and directional compartments (stick, zeppelin, tensor, NODDI [12]). The binary release also allows for the estimation of the DDI compartment model [9]. 

*Example:* this estimates a multi-compartment model with two components (-n 2), each anisotropic component being a tensor (-c 3), a free water compartment (-F) and an isotropic restricted water compartment (-R), using the variable projection method (--ml-mode 2) and the Levenberg Marquardt algorithm. 

.. prompt:: bash $

	animaMCMEstimator -i DWI.nii.gz -o MCM_n2.mcm -O MCM_n2_b0.nrrd -g DWI.bvec -b DWI.bval -n 2 -c 3 -F -R --optimizer levenberg --ml-mode 2

More options of this tool are available in the documentation of **animaMCMEstimator**. Interestingly the number of random restarts has a large influence on computation time (but also on precision). We fix it by default to 6 but getting to smaller values will make the estimation faster.

The output file format is the MCM format that may be read by `medInria <http://med.inria.fr>`_ from version 3.0. The format is an XML file linking to individual compartment images and weights.

*Example:*
This .mcm file references image files describing each compartment (together with a flag for its type). In addition, it references a vector image containing for each voxel the respective weights of each compartment.

::

	<?xml version="1.0"?>
	<Model>
	<Weights>testing_weights.nrrd</Weights>
	<Compartment>
	<Type>FreeWater</Type>
	<FileName>testing_0.nrrd</FileName>
	</Compartment>
	<Compartment>
	<Type>IRWater</Type>
	<FileName>testing_1.nrrd</FileName>
	</Compartment>
	<Compartment>
	<Type>Tensor</Type>
	<FileName>testing_2.nrrd</FileName>
	</Compartment>
	<Compartment>
	<Type>Tensor</Type>
	<FileName>testing_3.nrrd</FileName>
	</Compartment>

MCM model averaging
"""""""""""""""""""

In addition to estimation, we provide two ways of performing model selection and averaging:

* **animaMCMEstimator** proposes model selection based on AICc criterion (option -M): in that case, all models from 0 to N will be estimated and the one with the optimal AICc will be kept
* The binary Anima tools include one called **animaMCMModelAveraging** that implements the method proposed in [10]. This method uses outputs from model estimations from 0 to N fiber compartments and their AICc scores (produced by **animaMCMEstimator**) to compute an average MCM volume with simplification to the optimal number of fibres in each voxel. This tool is not yet open source and as such will be distributed only as binary versions in the Anima releases.

*Examples:*

* this computes a multi-tensor model at each voxel, with at most 3 anisotropic compartments per voxel (this number being decided based on the AICc criterion).

.. prompt:: bash $

	animaMCMEstimator -i DWI.nii.gz -o MCM_n3_MS.mcm -O MCM_n3_MS_b0.nrrd -g DWI.bvec -b DWI.bval -n 3 -c 3 -FR --optimizer levenberg --ml-mode 2 -M

* this performs model averaging as proposed in [10], with model simplification. It uses as an input two text files, each having on each line an image file name. *listMCM.txt* contains the list of MCM files to be averaged (from 0 to N anisotropic compartments) and *listAIC.txt* contains the corresponding AICc files (all these images may be written from **animaMCMEstimator**).

.. prompt:: bash $

	animaMCMModelAveraging -i listMCM.txt -o MCM_avg.mcm -a listAIC.txt -m MCM_avg_mose.nrrd -C

MCM processing
""""""""""""""

**animaMCMAverageImages** provides a way to average several volumes of MCM into just one (e.g. an atlas of those images), using the averaging and interpolation framework proposed in [11]. It works in a similar manner to the `animaAverageImages` described in the basic tools page.

DTI estimation and processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DTI estimation
""""""""""""""

DTI estimation is performed using two tools in ANIMA, implementing basic matrix-based DTI estimation and extrapolation.

**animaDTIEstimator** takes as inputs a 4D DWI image, a set of gradient directions and b-values and estimates tensors at each voxel. Gradient directions may be in the medInria format (one line per gradient) or the bvec format. B-values may be specified using a single number or a text file (either one line for each volume b-value or a bval file). Estimated tensors may be degenerated in some places. In that case, the tool outputs either zero values or the degenerated tensors depending on the -K option.

*Note:* In all Anima tools, the tensors are stored using a 6-component vector image representing the upper diagonal part of the tensors. These values are stored in column-first order.

DTI scalar maps
"""""""""""""""

**animaComputeDTIScalarMaps** computes the usual fractional anisotropy (FA), apparent diffusivity coefficient (ADC), axial (AD) and radial diffusivity (RD) maps from a tensor image.

Log-Euclidean tools
"""""""""""""""""""

These tools implement Arsigny et al. log and exponential maps on tensors:

* **animaLogTensors** computes the log map of tensors. The -S option switches between the vector representation and matrix representation of the log (sqrt(2) scaling factor on non diagonal terms).

* **animaExpTensors** computes the exponential map of log-vectors. The -S option is the equivalent of the one in **animaLogTensors**: it divides non diagonal values by sqrt(2).

ODF estimation and processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In all Anima tools, the ODFs are represented in the real spherical harmonics basis proposed by Descoteaux et al. in [2]. Coefficients are stored in vector images as explained in that publication.

ODF estimation
""""""""""""""

**animaODFEstimator** estimates ODFs at each voxel using one of two estimation methods: (1) Descoteaux et al. [2] with or without regularization, with or without ODF spherical deconvolution [3], and (2) Aganj et al. [4] providing naturally normalized ODFs at each voxel. The amount of ODF spherical deconvolution may be specified with the -s parameter, the estimation method with -R. 

*Example:* this estimates ODFs of order 6 from DWI.nii.gz using Aganj et al. method.

.. prompt:: bash $

	animaODFEstimator -i DWI.nii.gz -o ODF.nii.gz -g grads.bvec -k 6 -R

Generalized FA
""""""""""""""

**animaGeneralizedFA** computes the generalized fractional anisotropy from an image of ODFs stored in our format.

Tractography
------------

Anima implements tractography based on the three supported models: DTI, ODFs and MCM. It can be further divided into two classes of tractography methods: deterministic and probabilistic. All algorithms output fibers either in .vtk, .vtp (VTK format) or .fds (a meta-fibers format that can easily be read by `medInria <http://med.inria.fr>`_).

Deterministic tractography
^^^^^^^^^^^^^^^^^^^^^^^^^^

Deterministic tractography algorithms are described in [5]. They implement FACT [6] for DTI and a modified version of it to handle crossing fibers for ODFs. Those algorithms progress step by step following the local directions provided by the local model available and stopping if some criterions are met (local fiber angle, fiber length, FA threshold, ...).

* **animaDTITractography** implements DTI based deterministic tractography.
* **animaMCMTractography** implements multi-compartment models based deterministic tractography.

Probabilistic tractography
^^^^^^^^^^^^^^^^^^^^^^^^^^

Probabilistic tractography tools implement for MCM, ODF and DTI our multi-modal particle filtering framework for probabilistic tractography [7]. It relies on the simultaneous propagation of particles and their filtering relative to previous directions and the current model. This method further implements clustering of the particles to retain multi-modality, i.e. branching fibers.

* **animaDTIProbabilisticTractography** implements the filter for DTI tractography.
* **animaODFProbabilisticTractography** implements the filter for ODF tractography.
* **animaMCMProbabilisticTractography** implements the filter for multi-compartment models tractography.

Tractography tools
^^^^^^^^^^^^^^^^^^

**animaFibersApplyTransformSerie** works in the same way as resampler tools provided on the :doc:`registration page <registration>` except that it applies a series of transformations to a set of fibers. Please refer to that section for more details.

**animaFibersFilterer** uses a regions of interest (labeled) image to filter a set of fibers. The ROI image is a label image provided with the option -r. The -t and -f options can be given multiple times and are used to tell which labels a single fiber should go through (-t) and which labels should not be touched (-f).

*Example:* this filters the input fibers telling each fiber can be kept if it touches labels 1 and 2, but not 3.

.. prompt:: bash $

	animaFibersFilterer -i fibers.fds -o filtered_fibers.fds -r roi_image.nrrd -t 1 -t 2 -f 3 

References
----------

1. Vincent Arsigny, Pierre Fillard, Xavier Pennec, and Nicholas Ayache. *Log-Euclidean Metrics for Fast and Simple Calculus on Diffusion Tensors*. Magnetic Resonance in Medicine, 56(2):411-421, August 2006.
2. Descoteaux, M., Angelino, E., Fitzgibbons, S., Deriche, R. *Regularized, Fast, and Robust Analytical Q-Ball Imaging*. Magnetic Resonance in Medicine 58, 497–510, 2007.
3. Descoteaux M, Deriche R, Knösche TR, Anwander A. *Deterministic and probabilistic tractography based on complex fibre orientation distributions*. IEEE Transactions on Medical Imaging, 28(2):269-86, 2009.
4. Iman Aganj, Christophe Lenglet, Guillermo Sapiro, Essa Yacoub, Kamil Ugurbil, Noam Harel. *Reconstruction of the orientation distribution function in single‐and multiple‐shell q‐ball imaging within constant solid angle*. Magnetic Resonance in Medicine, 64(2):554-566, 2010.
5. Nicolas Wiest-Daesslé, Olivier Commowick, Aymeric Stamm, Patrick Perez, Christian Barillot, Romuald Seizeur, Sylvain Prima. *Comparison of 3 Diffusion Models to Track the Hand Motor Fibers within the Corticospinal Tract Using Functional, Anatomical and Diffusion MRI*. MICCAI 2011 Workshop on Computational Diffusion MRI (CDMRI'11), pp 150-157, Sep 2011.
6. Susumu Mori, Barbara J. Crain, V. P. Chacko, Peter C. M. Van Zijl. *Three-dimensional tracking of axonal projections in the brain by magnetic resonance imaging*. Annals of Neurology, 45(2):265–269, 1999.
7. Aymeric Stamm, Olivier Commowick, Christian Barillot, Patrick Perez. *Adaptive Multi-modal Particle Filtering for Probabilistic White Matter Tractography*. Information Processing in Medical Imaging, pp 594-606, 2013.
8. Aymeric Stamm, Olivier Commowick, Simon K. Warfield, Simone Vantini. *Comprehensive Maximum Likelihood Estimation of Diffusion Compartment Models Towards Reliable Mapping of Brain Microstructure*. 19th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2016.
9. Aymeric Stamm, Patrick Pérez, Christian Barillot. *A new multi-fiber model for low angular resolution diffusion MRI*. IEEE International Symposium on Biomedical Imaging, 2012.
10. Aymeric Stamm, Olivier Commowick, Patrick Pérez, Christian Barillot. *Fast Identification of Optimal Fascicle Configurations from Standard Clinical Diffusion MRI Using Akaike Information Criterion*. IEEE International Symposium on Biomedical Imaging, 2014.
11. Renaud Hédouin, Olivier Commowick, Aymeric Stamm, Christian Barillot. *Interpolation and Averaging of Multi-Compartment Model Images*, 18th International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 354-362, 2015.
12. Hui Zhang, Torben Schneider, Claudia A. Wheeler-Kingshott, Daniel C. Alexander. *NODDI: Practical in vivo neurite orientation dispersion and density imaging of the human brain*, NeuroImage, 61:4, 1000-1016, 2012.