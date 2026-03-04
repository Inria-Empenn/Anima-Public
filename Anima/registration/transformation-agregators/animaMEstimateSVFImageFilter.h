#pragma once

#include <iostream>
#include <itkImage.h>
#include <itkImageToImageFilter.h>

namespace anima {

template <class TScalarType, unsigned int NDegreesOfFreedom,
          unsigned int NDimensions = 3>
class MEstimateSVFImageFilter
    : public itk::ImageToImageFilter<
          itk::Image<itk::Vector<TScalarType, NDegreesOfFreedom>, NDimensions>,
          itk::Image<itk::Vector<TScalarType, NDegreesOfFreedom>,
                     NDimensions>> {
public:
  /** Standard class typedefs. */
  using Self = MEstimateSVFImageFilter;
  using WeightImageType = itk::Image<TScalarType, NDimensions>;
  using WeightImagePointer = typename WeightImageType::Pointer;
  using TInputImage =
      itk::Image<itk::Vector<TScalarType, NDegreesOfFreedom>, NDimensions>;
  using TOutputImage =
      itk::Image<itk::Vector<TScalarType, NDegreesOfFreedom>, NDimensions>;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MEstimateSVFImageFilter, ImageToImageFilter);

  using OutputPixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;
  using InputIndexType = typename TInputImage::IndexType;
  using InputPointType = typename TInputImage::PointType;

  /** Image typedef support */
  using InputImagePointer = typename TInputImage::Pointer;
  using OutputImagePointer = typename TOutputImage::Pointer;

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(WeightImage, WeightImagePointer);

  itkSetMacro(FluidSigma, double);
  itkSetMacro(MEstimateFactor, double);

  itkSetMacro(ConvergenceThreshold, double);
  itkSetMacro(MaxNumIterations, unsigned int);

protected:
  MEstimateSVFImageFilter() {
    m_FluidSigma = 4.0;
    m_MEstimateFactor = 1.0;
    m_AverageResidualValue = 1.0;
  }

  virtual ~MEstimateSVFImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MEstimateSVFImageFilter);

  bool checkConvergenceThreshold(OutputPixelType &outValOld,
                                 OutputPixelType &outVal);

  WeightImagePointer m_WeightImage;

  double m_FluidSigma, m_MEstimateFactor;
  std::vector<unsigned int> m_NeighborhoodHalfSizes;
  unsigned int m_MaxNumIterations;
  double m_ConvergenceThreshold;

  std::vector<double> m_InternalSpatialKernelWeights;
  std::vector<InputIndexType> m_InternalSpatialKernelIndexes;

  // Internal parameter
  double m_AverageResidualValue;
};

} // end namespace anima

#include "animaMEstimateSVFImageFilter.hxx"
