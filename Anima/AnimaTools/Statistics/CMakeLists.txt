cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME)
set(PROJECT_PREFIX anima${DIRNAME})

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(AverageImages)
add_subdirectory(BootStrap4dVolume)
add_subdirectory(BuildSamples)
add_subdirectory(CovarianceImages)
add_subdirectory(DirichletEstimation)
add_subdirectory(GammaEstimation)
add_subdirectory(LocalPatchCovarianceDistance)
add_subdirectory(LocalPatchMeanDistance)
add_subdirectory(LowMemoryTools)
add_subdirectory(RoiIntensitiesStats)
add_subdirectory(SampleImageFromDistribution)
add_subdirectory(StatisticalTests)
