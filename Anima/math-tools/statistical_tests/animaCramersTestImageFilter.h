#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <iostream>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <vector>

namespace anima {

template <class PixelScalarType>
class CramersTestImageFilter : public anima::MaskedImageToImageFilter<
                                   itk::VectorImage<PixelScalarType, 3>,
                                   itk::Image<PixelScalarType, 3>> {
public:
  /** Standard class typedefs. */
  using Self = CramersTestImageFilter<PixelScalarType>;
  using TInputImage = itk::VectorImage<PixelScalarType, 3>;
  using TOutputImage = itk::Image<PixelScalarType, 3>;
  using Superclass = anima::MaskedImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(CramersTestImageFilter, MaskedImageToImageFilter);

  /** Image typedef support */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using MaskImageType = typename Superclass::MaskImageType;
  using MaskImagePointer = typename MaskImageType::Pointer;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  /** Set the number of samples to build the distribution */
  itkSetMacro(NbSamples, unsigned long int);
  itkGetMacro(NbSamples, unsigned long int);

  void SetFirstGroupIndexes(std::vector<unsigned int> &group) {
    m_FirstGroup.clear();
    m_FirstGroupSize = group.size();

    for (unsigned int i = 0; i < m_FirstGroupSize; ++i)
      m_FirstGroup.push_back(group[i]);
  }

  void SetSecondGroupIndexes(std::vector<unsigned int> &group) {
    m_SecondGroup.clear();
    m_SecondGroupSize = group.size();

    for (unsigned int i = 0; i < m_SecondGroupSize; ++i)
      m_SecondGroup.push_back(group[i]);
  }

  void AddOutlierMask(MaskImageType *mask) { m_OutlierMasks.push_back(mask); }

protected:
  CramersTestImageFilter() : Superclass() {
    m_NbSamples = 5000;
    m_FirstGroup.clear();
    m_SecondGroup.clear();
    m_FirstGroupSize = 0;
    m_SecondGroupSize = 0;

    m_UseOutlierMasks = false;
    m_OutlierMasks.clear();
  }

  virtual ~CramersTestImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

  //! Bootstrap samples generator for group comparison
  void GenerateBootStrapSamples();

  //! Actually bootstraps a p-value from the computed distance matrix and group
  //! samples
  double BootStrap(vnl_matrix<double> &groupDistMatrix,
                   std::vector<double> &inlierWeights);

  //! Computes the Cramers' statistic knowing the distance matrix and the two
  //! groups
  double CramerStatistic(vnl_matrix<double> &grpDistMatrix,
                         std::vector<double> &inlierWeights,
                         unsigned int index);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(CramersTestImageFilter);

  std::vector<unsigned int> m_FirstGroup, m_SecondGroup;
  std::vector<std::vector<unsigned int>> m_SamplesFirstGroup,
      m_SamplesSecondGroup;

  unsigned long int m_NbSamples;
  unsigned int m_FirstGroupSize, m_SecondGroupSize;

  std::vector<MaskImagePointer> m_OutlierMasks;
  bool m_UseOutlierMasks;
};

} // end namespace anima

#include "animaCramersTestImageFilter.hxx"
