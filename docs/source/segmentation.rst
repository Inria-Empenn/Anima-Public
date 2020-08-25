Segmentation tools
==================

This part describes segmentation features in Anima. It, as of now, includes multiple sclerosis segmentation tools and graph-cut segmentation. 

Basic segmentation tools
------------------------

We provide a set of basic tools related to segmentation that can be used easily.

Dice score computation
^^^^^^^^^^^^^^^^^^^^^^

**animaDiceMeasure** computes the Dice score between two labeled segmentations. For binary images, it computes the Dice or Jaccard scores. For multi-labeled images, it either outputs the individual label scores or the total overlap score as in Klein et al. [1] (``-T`` option).

Image thresholding
^^^^^^^^^^^^^^^^^^

**animaThrImage** thresholds an image, using one or several thresholds. 

*Example:* this thresholds the input Image.nii.gz to keep only voxels whose values are between 0.5 and 1.5.

.. code-block:: sh

	animaThrImage -i Image.nii.gz -t 0.5 -u 1.5 -o BinaryImage.nii.gz

**animaOtsuThrImage** computes the Otsu automatic thresholding in several classes of an input image.

Image holes filler
^^^^^^^^^^^^^^^^^^

**animaFillHoleImage** provides a simple tool that fills small holes in segmentations that are not connected with the background.

Image masking
^^^^^^^^^^^^^

**animaMaskImage** simply takes an image and a binary mask and applies the mask to the input image.

Isosurface extraction
^^^^^^^^^^^^^^^^^^^^^

**animaIsosurface** takes as an input image and extracts its isosurface at the level given by the ``-t`` option. Decimation and smoothing of the resulting surface can also be controlled.

Graph cut segmentation
----------------------

We provide an implementation of Lecoeur et al. graph cut segmentation algorithm [1], using the spectral gradient for multimodal images. It is implemented in the **animaGraphCut** tool.

Multiple sclerosis lesions segmentation
---------------------------------------

We provide an implementation of one of the multiple sclerosis (MS) segmentation algorithms developed in the Visages team [2,3]. This method is based on a tissue classification of multiple channel images and uses graph cut to delineate T2 hyper-intense lesions. It takes as an input three modalities among T1, T2, DP and FLAIR as well as a brain mask.

*Example:* this command line uses T1, T2 and FLAIR images as an input, as well as brain_mask.nii.gz as a skull-stripping mask and compute their GCEM segmentation.

.. code-block:: sh

	animaGcStremMsLesionsSegmentation -m brain_mask.nii.gz -a 1 --rej 0.2 --min-fuzzy 1 --max-fuzzy 2 --intT2 3 --intFLAIR 2 --rb --ml 3 -i T1.nii.gz -j T2.nii.gz -l FLAIR.nii.gz --ini 2 -o lesion_seg.nii.gz --out-csf csf_seg.nii.gz --out-gm gm_seg.nii.gz --out-wm wm_seg.nii.gz --out-gc gc_seg.nii.gz 

The parameters are as follows: 

* ``-a`` option specifies the type of algorithm (0 is STREM [2], 1 is GCEM [3]). 
* ``--rej`` : percentage of outliers rejection
* ``--min-fuzzy`` : minimal value for fuzzy rules
* ``--max-fuzzy`` : maximal value for fuzzy rules
* ``--ml`` : minimal lesion size (mm3)
* ``--ini`` : initialization method (0: atlas, 1: based on DP, 2: based on FLAIR)
* ``-o`` : output lesion segmentation
* ``--out-csf`` : output CSF map
* ``--out-gm`` : output GM map
* ``--out-wm`` : output WM map

Segmentation validation tools
-----------------------------

We now provide tools that were used for the validation of the `MS segmentation challenge <http://go.nature.com/2SW1DhA>`_ held in 2016 and that can be used for other tasks as well.

Dice measure
^^^^^^^^^^^^

Anima comes with two tools for computing the Dice measure:

* **animaDiceMeasure** computes Dice scores between two label images (i.e. with each pixel having an integer label). It can either output the Dice scores for each label individually or compute the total overlap score as proposed by Klein et al. [5] (``-T`` option)
* **animaFuzzyDiceMeasure** computes the generalized Dice score between fuzzy segmentations of a structure as proposed by Crum et al. [6]

Both tools have an option to output the Jaccard score instead of the Dice score (``-J`` option), both scores being related by a monotonic function.

Detected components
^^^^^^^^^^^^^^^^^^^

**animaDetectedComponents** provides an evaluation tool that counts the number of connected components in a reference binary segmentation that are actually detected by a test segmentation. This does not mean that the components are very well segmented but rather that the test segmentation detect sufficiently the components. The output is a CSV file detailing for the reference image the number of connected components, their sizes and if they were detected by the test segmentation.

Segmentation performance analyzer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**animaSegPerfAnalyzer** includes all metrics that were computed for the MS segmentation challenge at MICCAI in 2016. This tool was used to compute the results available as supplementary material `here <https://doi.org/10.5281/zenodo.1307652>`_. Please refer to our article [4] for more details on the measures.

References
----------

1. Jérémy Lecoeur, Sean Patrick Morrissey, Jean-Christophe Ferré, Douglas Arnold, D. Louis Collins, Christian Barillot. *Multiple Sclerosis Lesions Segmentation using Spectral Gradient and Graph Cuts*. Medical Image Analysis on Multiple Sclerosis (validation and methodological issues), MICCAI workshop, 2008.
2. Daniel García-Lorenzo, Sylvain Prima, Douglas Arnold, Louis Collins, Christian Barillot. *Trimmed-likelihood estimation for focal lesions and tissue segmentation in multisequence MRI for multiple sclerosis*. IEEE Transactions on Medical Imaging, 30 (8), pp.1455-67, 2011.
3. Daniel García-Lorenzo, Jérémy Lecoeur, Douglas Arnold, D. Louis Collins, Christian Barillot. *Multiple Sclerosis lesion segmentation using an automatic multimodal Graph Cuts*. 12th International Conference on Medical Image Computing and Computer Assisted Intervention, LNCS 5762, pp.584-591, 2009.
4. O\. Commowick et al\. *Objective Evaluation of Multiple Sclerosis Lesion Segmentation using a Data Management and Processing Infrastructure*. Scientific Reports, 8(1), 2018
5. Klein, A, Andersson, J, Ardekani, BA, Ashburner, J, Avants, B, Chiang, M-C, Christensen, GE, Collins, DL, Gee, J, Hellier, P, Song, JH, Jenkinson, M, Lepage, C, Rueckert, D, Thompson, P, Vercauteren, T, Woods, RP, Mann, JJ, Parsey, RV. *Evaluation of 14 nonlinear deformation algorithms applied to human brain MRI registration*. NeuroImage. 46(3): 786-802. 2009.
6.  W.R. Crum, O. Camara and D.L.G. Hill. *Generalized Overlap Measures for Evaluation and Validation in Medical Image Analysis*. IEEE Transactions on Medical Imaging. 25(11):1451-1461. 2006.

