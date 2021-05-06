Relaxometry tools
=================

Estimation
----------

T1 relaxation
^^^^^^^^^^^^^

We include an implementation of the DESPOT1 [1] algorithm for T1 relaxation time estimation: **animaT1RelaxometryEstimation**. It optionally can take a B1 inhomogeneity map as an input to correct for it (``-b`` parameter). It outputs two images: T1 relaxation times (``-o`` output) and M0 map (proportional to proton density, ``-O`` output).

T2 relaxation
^^^^^^^^^^^^^

We include two ways to estimate T2 relaxation times from T2 relaxometry sequences. Those two algorithms suppose acquisitions were made with an equal echo spacing that is specified to the algorithm. Those two tools are:

* **animaT2RelaxometryEstimation** estimates T2 relaxation using the regular exponential decay equation with a single T2 value. It may take an input T1 map for better estimation, and outputs M0 and T2 maps.
* **animaT2EPGRelaxometryEstimation** estimates T2 relaxation using the EPG algorithm with a single T2 value, to account for stimulated echoes due to B1 inhomogeneity [2]. It may take input T1 and B1 maps for better estimation, and outputs M0 and T2 maps.

Multi-compartment T2 estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide from Anima 3.0 new tools for multi-compartment T2 and myelin water fraction estimation. We also provide implementations of previous methods from the literature on multi-peak T2 estimation. While this second approach is appealing, we find it more stable to use multi-compartment T2 estimation rather than multi-peak and rather provide those for comparison purposes.

* **animaGammaMixtureT2RelaxometryEstimation** implements variable projection estimation of the parameters and weights of three T2 Gamma distributions using two different modes (toggled by the ``-U`` option) [3,4]: only middle T2 compartment mean estimation or all class mean parameters estimation. In both cases, all weights are deduced from the estimated parameters using variable projection
* **animaGMMT2RelaxometryEstimation** implements, for robustness to clinical acquisitions, a fixed parameter estimation of a Gaussian T2 mixture [7]. In this implementation, only the weights of the three T2 compartments are estimated, their PDFs being fixed according to prior knowledge on the tissues.
* **animaMultiT2RelaxometryEstimation** provides an implementation of several methods of the literature for multi-peak T2 estimation [2,5,6]. It provides several types of regularization: Tikhonov, Laplacian or non local regularization. Again these methods make the deduction of myelin water fraction more difficult and quite sensitive to the regularization.

MRI simulation
--------------

Several MR simulation tools are included in ANIMA, which simulate sequences from relaxation time maps. All of them are described in detail `here <https://team.inria.fr/empenn/files/2017/08/mr_simulation_guide.pdf>`_. 

References
----------

1. Deoni, S. *High-resolution T1 mapping of the brain at 3T with driven equilibrium single pulse observation of T1 with high-speed incorporation of RF field inhomogeneities (DESPOT1-HIFI)*. J Magn Reson Imaging 26(4):1106–1111, 2007.
2. Kelvin J Layton, Mark Morelande, David Wright, Peter M Farrell, Bill Moran, Leigh Johnston. *Modelling and estimation of multicomponent T2 distributions*. IEEE transactions on medical imaging, 32(8):1423–34, 2013.
3. Sudhanya Chatterjee, Olivier Commowick, Onur Afacan, Simon K. Warfield, Christian Barillot. *Multi-Compartment Model of Brain Tissues from T2 Relaxometry MRI Using Gamma Distribution*. ISBI 2018.
4. Sudhanya Chatterjee, Olivier Commowick, Simon K. Warfield, Christian Barillot. *Multi-Compartment T2 Relaxometry Model Using Gamma Distribution Representations: A Framework for Quantitative Estimation of Brain Tissue Microstructures*. ISMRM 2017.
5. Prasloski et al. *Applications of stimulated echo correction to multicomponent T2 analysis*. MRM, 67(6):1803-1814, 2012.
6. Yoo et al. *Non-local spatial regularization of MRI T2 relaxation images for myelin water quantification*. MICCAI, pp 614-621, 2013.
7. S\. Chatterjee, O\. Commowick, O\. Afacan, B\. Combes, A\. Kerbrat, S\.K\. Warfield, C\. Barillot. *A 3-year follow-up study of enhancing and non enhancing multiple sclerosis lesions in MS patients with clinically isolated syndrom using a multi-compartment T2 relaxometry model*, ISMRM, 2018.
