#pragma once

#include <iostream>
#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkObject.h>
#include <itkVector.h>
#include <itkVectorImage.h>

namespace anima {

template <class TInputImage>
class NonLocalMeansTemporalImageFilter
    : public itk::ImageToImageFilter<TInputImage, TInputImage> {
public:
  /** Define pixel types  */
  using InputPixelType = typename TInputImage::PixelType;
  using OutputPixelType = double;

  /** Convenient typedefs for simplifying declarations. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputImageIndexType = typename InputImageType::IndexType;

  using OutputImageType = InputImageType;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputImageRegionType = typename OutputImageType::RegionType;

  /** Standard "Self" & Superclass typedef. */
  using Self = NonLocalMeansTemporalImageFilter;
  using Superclass = itk::ImageToImageFilter<InputImageType, OutputImageType>;

  /** SmartPointer typedef support  */
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory.  */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NonLocalMeansTemporalImageFilter, itk::ImageToImageFilter);

  /** Extract dimension from input image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Noise used to determine weightes */
  enum WEIGHT { EXP, RICIAN };

  itkSetMacro(PatchHalfSize, unsigned int);
  itkSetMacro(SearchNeighborhood, unsigned int);
  itkSetMacro(SearchStepSize, unsigned int);
  itkSetMacro(WeightThreshold, double);
  itkSetMacro(BetaParameter, double);
  itkSetMacro(MeanMinThreshold, double);
  itkSetMacro(VarMinThreshold, double);
  itkSetMacro(WeightMethod, WEIGHT);

protected:
  NonLocalMeansTemporalImageFilter()
      : m_MeanMinThreshold(0.95), m_VarMinThreshold(0.5),
        m_WeightThreshold(0.0), m_BetaParameter(1.0), m_PatchHalfSize(1),
        m_SearchStepSize(1), m_SearchNeighborhood(5), m_WeightMethod(EXP) {}

  virtual ~NonLocalMeansTemporalImageFilter() {}

  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(NonLocalMeansTemporalImageFilter);

  double m_MeanMinThreshold;
  double m_VarMinThreshold;
  double m_WeightThreshold;
  double m_BetaParameter;

  unsigned int m_PatchHalfSize;
  unsigned int m_SearchStepSize;
  unsigned int m_SearchNeighborhood;
  WEIGHT m_WeightMethod;
};

} // end of namespace anima

#include "animaNonLocalMeansTemporalImageFilter.hxx"
