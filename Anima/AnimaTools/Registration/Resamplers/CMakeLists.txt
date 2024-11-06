add_subdirectory(ApplyTransformSerie)
add_subdirectory(ApplyVtkTransform)
add_subdirectory(ApplyDistortionCorrection)

if (USE_VTK AND VTK_FOUND)
  add_subdirectory(FibersApplyTransformSerie)
endif()

add_subdirectory(MCMApplyTransformSerie)
add_subdirectory(ODFApplyTransformSerie)
add_subdirectory(TensorApplyTransformSerie)
add_subdirectory(PyramidImage)
