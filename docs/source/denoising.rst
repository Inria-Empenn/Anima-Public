Denoising tools
===============

Noise generation
----------------

**animaNoiseGenerator** adds Gaussian or Rician noise to a scalar image at a user specified SNR (an external reference input image is necessary to compute the average reference signal value). The -n may be used to output multiple samples of noisy images from the input.

*Example:* this adds Gaussian noise (relative to the average signal computed from refSignalImage.nii.gz) to Image.nii.gz. It will provide as outputs 10 images named OutputImage\_??.nii.gz

.. prompt:: bash $

	animaNoiseGenerator -i Image.nii.gz -o OutImage.nii.gz -G -n 10 -r refSignalImage.nii.gz

NL-Means denoising
------------------

NL-Means denoising is an implementation of the two papers from the team to denoise images corrupted by Gaussian [1] or Rician noise [2].

NL-Means for scalar images
^^^^^^^^^^^^^^^^^^^^^^^^^^

**animaNLMeans** implements patch-based denoising for scalar images, considering it as a single image without a temporal dimension. It includes all parameters described in [1], including patch half size (real size being 2N+1), patch half search neighborhood, beta parameter for weighting, ... One can choose between [1] and [2] with the -W option.

*Example:* this performs denoising of Image.nii.gz with the default parameters, a beta parameter of 0.5 and assuming Rician noise.

.. prompt:: bash $

	animaNLMeans -i Image.nii.gz -o OutImage.nii.gz -W 1 -b 0.5

NL-Means for images with a temporal dimension
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**animaNLMeansTemporal** works in the same way as **animaNLMeans** except that it considers the last dimension of the input image as a temporal dimension and therefore processes each temporal volume independently (useful e.g. for DWI images).

References
----------

1. Pierrick Coupé, Pierre Yger, Sylvain Prima, Pierre Hellier, Charles Kervrann, Christian Barillot. *An optimized blockwise nonlocal means denoising filter for 3-D magnetic resonance images*. IEEE Transactions on Medical Imaging. 27(4): 425-441. 2008.
2. Nicolas Wiest-Daesslé, Sylvain Prima, Pierrick Coupé, Sean Patrick Morrissey, Christian Barillot. *Rician noise removal by non-Local Means filtering for low signal-to-noise ratio MRI: applications to DT-MRI*. 11th International Conference on Medical Image Computing and Computer-Assisted Intervention. Springer, 5242 (Pt 2), pp.171-179, 2008.