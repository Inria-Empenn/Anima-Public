
## ##############################################################################
## Here go the add_subdirectories, no code should be at the root
## ##############################################################################

add_subdirectory(Common)
add_subdirectory(Diffusion)
add_subdirectory(Filtering)
add_subdirectory(Interpolator)
add_subdirectory(Registration)
add_subdirectory(Regularization)
add_subdirectory(Relaxometry)
add_subdirectory(Segmentation)



add_subdirectory(Common/NonLocalBase)

add_subdirectory(Diffusion/dti/BackgroundNoiseVarianceEstimator)
add_subdirectory(Diffusion/LogTensorImageFilter)
add_subdirectory(Diffusion/MCM)
add_subdirectory(Diffusion/MCMSimilarity)
add_subdirectory(Diffusion/Tractography)

add_subdirectory(Filtering/dti/filters)
add_subdirectory(Filtering/ClusteringFilter)
add_subdirectory(Filtering/Jacobian)

#add_subdirectory(Interpolator)

add_subdirectory(Registration/Algorithms)
add_subdirectory(Registration/Common)
add_subdirectory(Registration/DiffusionMetrics)
add_subdirectory(Registration/GenericMetrics)
add_subdirectory(Registration/PyramidalBMRegistration)
add_subdirectory(Registration/PyramidalDenseSVFMatchingBridge)
add_subdirectory(Registration/ResamplersFilters)
add_subdirectory(Registration/SVFOperations)
add_subdirectory(Registration/TransformationAgregators)
add_subdirectory(Registration/Transformations)

add_subdirectory(Regularization/Filters)

add_subdirectory(Relaxometry/Estimation)
add_subdirectory(Relaxometry/Simulation)

# add_subdirectory(Segmentation)

