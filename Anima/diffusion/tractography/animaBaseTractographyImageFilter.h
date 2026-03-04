#pragma once

#include <itkImage.h>
#include <itkVectorImage.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkProcessObject.h>
#include <mutex>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "AnimaTractographyExport.h"

#include <vector>

namespace anima {

class ANIMATRACTOGRAPHY_EXPORT BaseTractographyImageFilter
    : public itk::ProcessObject {
public:
  /** SmartPointer typedef support  */
  using Self = BaseTractographyImageFilter;
  using Superclass = itk::ProcessObject;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkTypeMacro(BaseTractographyImageFilter, itk::ProcessObject);

  using ModelImageType = itk::VectorImage<double, 3>;
  using ModelImagePointer = ModelImageType::Pointer;
  using VectorType = ModelImageType::PixelType;
  using RegionType = ModelImageType::RegionType;
  using ContinuousIndexType = itk::ContinuousIndex<double, 3>;

  using MaskImageType = itk::Image<unsigned short, 3>;
  using MaskImagePointer = MaskImageType::Pointer;
  using PointType = MaskImageType::PointType;
  using IndexType = MaskImageType::IndexType;

  using FiberType = std::vector<PointType>;
  using FiberVectorType = std::vector<FiberType>;

  typedef struct {
    BaseTractographyImageFilter *trackerPtr;
    std::vector<FiberVectorType> resultFibersFromThreads;
  } trackerArguments;

  virtual void SetInputImage(ModelImageType *input) { m_InputImage = input; }
  ModelImageType *GetInputImage() { return m_InputImage; }

  void SetSeedingMask(MaskImageType *mask) { m_SeedingImage = mask; }
  void SetFilteringMask(MaskImageType *mask) { m_FilteringImage = mask; }
  void SetForbiddenMask(MaskImageType *mask) { m_ForbiddenMaskImage = mask; }
  void SetCutMask(MaskImageType *mask) { m_CutMaskImage = mask; }

  void SetNumberOfFibersPerPixel(unsigned int num) {
    m_NumberOfFibersPerPixel = num;
  }

  void SetStepProgression(double num) { m_StepProgression = num; }

  void SetMaxFiberAngle(double num) { m_MaxFiberAngle = num; }
  void SetMinLengthFiber(double num) { m_MinLengthFiber = num; }
  void SetMaxLengthFiber(double num) { m_MaxLengthFiber = num; }

  void Update() ITK_OVERRIDE;

  void SetComputeLocalColors(bool flag) { m_ComputeLocalColors = flag; }
  void createVTKOutput(std::vector<std::vector<PointType>> &filteredFibers);
  vtkPolyData *GetOutput() { return m_Output; }

protected:
  BaseTractographyImageFilter();
  virtual ~BaseTractographyImageFilter() {}

  static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreadTracker(void *arg);
  void ThreadTrack(unsigned int numThread,
                   std::vector<FiberType> &resultFibers);
  void ThreadedTrackComputer(unsigned int numThread,
                             std::vector<FiberType> &resultFibers,
                             unsigned int startSeedIndex,
                             unsigned int endSeedIndex);

  void ComputeFiber(FiberType &fiber, PointType &initialDirection,
                    itk::ThreadIdType threadId);

  virtual void PrepareTractography();
  std::vector<FiberType> FilterOutputFibers(std::vector<FiberType> &fibers);

  virtual bool CheckModelCompatibility(VectorType &modelValue,
                                       itk::ThreadIdType threadId) = 0;

  using ImageBaseType = itk::ImageBase<ModelImageType::ImageDimension>;
  bool CheckIndexInImageBounds(IndexType &index, ImageBaseType *testImage);
  virtual bool CheckIndexInImageBounds(ContinuousIndexType &index) = 0;

  //! Computes value of model from data. May use SNR and previous model value to
  //! perform smart interpolation. Replaces SNR and modelValue by the outputs
  virtual void GetModelValue(ContinuousIndexType &index,
                             VectorType &modelValue) = 0;

  virtual std::vector<PointType>
  GetModelPrincipalDirections(VectorType &modelValue, bool is2d,
                              itk::ThreadIdType threadId) = 0;
  virtual PointType GetNextDirection(PointType &previousDirection,
                                     VectorType &modelValue, bool is2d,
                                     itk::ThreadIdType threadId) = 0;

  //! Computes new fiber point using Runge Kutta integration (better spread of
  //! fibers than Euler integration)
  virtual void ComputeNewFiberPoint(PointType &oldPoint,
                                    PointType &newDirection,
                                    PointType &newPoint,
                                    itk::ThreadIdType threadId);

  virtual void ComputeAdditionalScalarMaps() {}
  bool isZero(VectorType &value);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(BaseTractographyImageFilter);

  unsigned int m_NumberOfFibersPerPixel;

  double m_StepProgression;
  double m_MaxFiberAngle;

  double m_MinLengthFiber;
  double m_MaxLengthFiber;

  ModelImagePointer m_InputImage;
  MaskImagePointer m_SeedingImage, m_FilteringImage, m_ForbiddenMaskImage,
      m_CutMaskImage;

  std::vector<ContinuousIndexType> m_PointsToProcess;
  std::vector<unsigned int> m_FilteringValues;

  bool m_ComputeLocalColors;
  vtkSmartPointer<vtkPolyData> m_Output;

  std::mutex m_LockHighestProcessedSeed, m_LockProcessedPoints;
  unsigned int m_NumberOfProcessedPoints;
  int m_HighestProcessedSeed;
};

} // end of namespace anima
