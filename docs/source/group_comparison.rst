Group comparison tools
======================

Patient to group comparison
---------------------------

General patient to group comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Patient to group comparison tools provides ways to compare a single patient image (scalar, log-tensor fields, vector fields, ODFs,...) to a group (population) of control images. The most general executable for it is **animaPatientToGroupComparison**. It handles any vector image type that represents data in a vector space (e.g. scalar images, log-tensor images obtained from **animaLogTensors**, stacked multimodal data...). It implements the test proposed in [1] but without direction sampling of the data. It also handles dimensionality reduction through PCA on the control database (``-E`` and ``-e`` options).

*Example:* this tests at each voxel testedImage against the list of controls (specified in listControls.txt: one line = one image to load). For speed reasons, it is highly desirable to use computationMask.nii.gz, which can be the intersection of brain masks of the controls and the tested image. It outputs p-values and z-scores of the test.

.. code-block:: sh

	animaPatientToGroupComparison -i testedImage.nii.gz -I listControls.txt -o zScore.nii.gz -O pValues.nii.gz -m computationMask.nii.gz

ODF patient to group comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**animaPatientToGroupODFComparison** implements the same test as in [1] but for ODF data. In addition to parameters of **animaPatientToGroupComparison**, the user has to specify a set of gradient directions (bvec or text file) on which the ODF will be sampled.

Non local patient to group comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Non-local patient to group comparison tools implement the method proposed in [2] for DTI images. These consist of three tools. The main one is **animaNLMeansPatientToGroupComparison** which performs the comparison itself by looking for patches around a given position, similarly to NL-Means denoising. It adds four options (``-d``, ``-D``, ``--dbcovave``, ``--dbcovstd``) to perform the average and standard deviation tests on the test image. Those four options require images as an input. They are provided by two other tools: 

* **animaLocalPatchMeanDistance** takes as an input the list of controls and provides two outputs: the one from ``-o`` is to be used for ``-d`` input of **animaNLMeansPatientToGroupComparison**, the one from ``-O`` is to be used for ``-D`` input of **animaNLMeansPatientToGroupComparison**)
* **animaLocalPatchCovarianceDistance** takes as an input the list of controls and provides two outputs: the one from ``-o`` is to be used for ``--dbcovave`` input of **animaNLMeansPatientToGroupComparison**, the one from ``-O`` is to be used for ``--dbcovstd`` input of **animaNLMeansPatientToGroupComparison**)

Group comparison tools
----------------------

We provide an implementation of a population comparison tool proposed by Whitcher et al. [3] and used for neonates comparison in [4]. This test is based on the computation of inter- and intra-group distances computation and their use to derive a statistical test of groupies voxel differences. As such our implementation uses as an input vector images for which a L2 distance may be computed (such as log-vectors of tensors or scalar images) and outputs a map of p-values (not corrected for multiple comparisons (see next section).

*Example:* this computes the voxelwise Cramers tests p-values given its inputs. dataList.txt lists on each line the input files of all groups, firstGroup.txt and secondGroup.txt lists the indexes of patients in each group, starting with index 0.

.. code-block:: sh

	animaCramersTest -l dataList.txt -m computationMask.nrrd -o outputPValues.nrrd -f firstGroup.txt -s secondGroup.txt

Multiple comparisons correction
-------------------------------

When doing the above mentioned tests, multiple comparisons are being made that need to be corrected for. **animaFDRCorrectPValues** implements FDR correction as presented by Benjamini and Hochberg [5]. It provides as an output thresholded p-values at q-value specified by the ``-q`` option.

Low memory tools
----------------

All preceding tools have their **animaLowMem** counterparts, as they usually require a lot of memory to run which may not be available on a single computer. These low memory tools include options to split the input images in every direction and process either a single sub-part of them or all of them in a sequential order. 

**animaMergeBlockImages** can then use the output description file and a geometry image to rebuild the final result. 

*Warning:* these low memory tools are much slower as they spend most of their time reading images. However, they can be run on a cluster if available.

References
----------

1. Olivier Commowick, Adil Maarouf, Jean-Christophe Ferré, Jean-Philippe Ranjeva, Gilles Edan, Christian Barillot. *Diffusion MRI abnormalities detection with orientation distribution functions: A multiple sclerosis longitudinal study*. Medical Image Analysis, 22(1):114-123, 2015.
2. Olivier Commowick, Aymeric Stamm. *Non-local robust detection of DTI white matter differences with small databases*. 15th International Conference on Medical Image Computing and Computer Assisted Intervention, Oct 2012, Nice, France. 15 (Pt 3), pp.476-84, 2012, LNCS.
3. B\. Whitcher, J\.J\. Wisco, N\. Hadjikhani and D\.S\. Tuch\. *Statistical group comparison of diffusion tensors via multivariate hypothesis testing*. Magnetic Resonance in Medicine, 57(6):1065-1074, June 2007.
4. O\. Commowick, N\. I\. Weisenfeld, H\. Als, G\. B\. McAnulty, S\. Butler, L\. Lightbody, R\. M\. Robertson and S\. K\. Warfield. *Evaluation of White Matter in Preterm Infants With Fetal Growth Restriction*, In Proceedings of the Workshop on Image Analysis for the Developing Brain, held in conjunction with MICCAI'09, September 2009.
5. Y\. Benjamini and Y\. Hochberg. *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 57, No. 1 (1995), pp. 289-300.