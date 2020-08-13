Relaxometry scripts
===================

This section describes Anima scripts for relaxometry processing. 

T2 relaxometry script
----------------------------------

This script combines several processing steps to compute from a 4D T2 relaxomaetry sequence (plus possible additional images) T2 maps or multi-T2 maps following [1]. This script performs the following preprocessing:

* brain masking either from the provided T1-weighted high resolution image provided using the **-m** option or by default the first volume of the input 4D image
* computation of the mono T2 maps with B1 correction using **animaT2EPGRelaxometryEstimation** [2]
* computation of multi-T2 weight maps (modeled as three Gaussians as in [1]) as well as the myelin water fraction using **animaGMMT2RelaxometryEstimation** [1]

Some of these steps may be discarded using the **--no-\*** options available in the script help. If provided, the EPG algorithms may use an external T1 image provided by the user with the **-t** option. Results of this script are specified using the **-o** and **-g** options.

References
----------

1. S. Chatterjee, O. Commowick, O. Afacan, B. Combes, A. Kerbrat, S.K. Warfield, C. Barillot. *A 3-year follow-up study of enhancing and non enhancing multiple sclerosis lesions in MS patients with clinically isolated syndrom using a multi-compartment T2 relaxometry model*, ISMRM, 2018.
2. Kelvin J Layton, Mark Morelande, David Wright, Peter M Farrell, Bill Moran, Leigh Johnston. *Modelling and estimation of multicomponent T2 distributions*. IEEE transactions on medical imaging, 32(8):1423–34, 2013.
