add_subdirectory(ApplyDistortionCorrection)
add_subdirectory(ApplyTransformSerie)
add_subdirectory(ApplyVtkTransform)

if (USE_VTK AND VTK_FOUND)
  add_subdirectory(FibersApplyTransformSerie)
endif()

add_subdirectory(MCMApplyTransformSerie)
add_subdirectory(ODFApplyTransformSerie)
add_subdirectory(PyramidImage)
add_subdirectory(TensorApplyTransformSerie)
