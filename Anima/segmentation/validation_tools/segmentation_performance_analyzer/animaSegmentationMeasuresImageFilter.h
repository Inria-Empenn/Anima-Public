#pragma once

#include <animaNumberedThreadImageToImageFilter.h>
#include <itkNumericTraits.h>

namespace anima {

/** \class LabelSetMeasures
 * \brief Metrics stored per label
 * \ingroup ITKImageStatistics
 */
class SegPerfLabelSetMeasures {
public:
  // default constructor
  SegPerfLabelSetMeasures() {
    m_Source = 0;
    m_Target = 0;
    m_Union = 0;
    m_TrueNegative = 0;
    m_Intersection = 0;
    m_SourceComplement = 0;
    m_TargetComplement = 0;
  }

  // added for completeness
  SegPerfLabelSetMeasures &operator=(const SegPerfLabelSetMeasures &l) {
    if (this != &l) {
      m_Source = l.m_Source;
      m_Target = l.m_Target;
      m_Union = l.m_Union;
      m_TrueNegative = l.m_TrueNegative;
      m_Intersection = l.m_Intersection;
      m_SourceComplement = l.m_SourceComplement;
      m_TargetComplement = l.m_TargetComplement;
    }
    return *this;
  }

  unsigned long m_Source;
  unsigned long m_Target;
  unsigned long m_Union;
  unsigned long m_TrueNegative;
  unsigned long m_Intersection;
  unsigned long m_SourceComplement;
  unsigned long m_TargetComplement;
};

template <typename TLabelImage>
class SegmentationMeasuresImageFilter
    : public anima::NumberedThreadImageToImageFilter<TLabelImage, TLabelImage> {
public:
  /** Standard Self typedef */
  using Self = SegmentationMeasuresImageFilter;
  using Superclass =
      anima::NumberedThreadImageToImageFilter<TLabelImage, TLabelImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(SegmentationMeasuresImageFilter,
               anima::NumberedThreadImageToImageFilter);

  /** Image related typedefs. */
  using LabelImageType = TLabelImage;
  using LabelImagePointer = typename TLabelImage::Pointer;
  using LabelImageConstPointer = typename TLabelImage::ConstPointer;

  using RegionType = typename TLabelImage::RegionType;
  using SizeType = typename TLabelImage::SizeType;
  using IndexType = typename TLabelImage::IndexType;

  using LabelType = typename TLabelImage::PixelType;

  /** Type to use for computations. */
  using RealType = typename itk::NumericTraits<LabelType>::RealType;

  /** Type of the map used to store data per label */
  using MapType = std::map<LabelType, SegPerfLabelSetMeasures>;
  using MapIterator = typename MapType::iterator;
  using MapConstIterator = typename MapType::const_iterator;

  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TLabelImage::ImageDimension);

  /** Set the source image. */
  void SetSourceImage(const LabelImageType *image) {
    this->SetNthInput(0, const_cast<LabelImageType *>(image));
  }

  /** Set the target image. */
  void SetTargetImage(const LabelImageType *image) {
    this->SetNthInput(1, const_cast<LabelImageType *>(image));
  }

  /** Get the source image. */
  const LabelImageType *GetSourceImage(void) { return this->GetInput(0); }

  /** Get the target image. */
  const LabelImageType *GetTargetImage(void) { return this->GetInput(1); }

  /** Get the label set measures */
  MapType GetLabelSetMeasures() { return this->m_LabelSetMeasures; }

  /**
   * tric overlap measures
   */
  /** measures over all labels */
  RealType getSensitivity();
  RealType getSpecificity();
  RealType getPPV();
  RealType getNPV();
  RealType getUnionOverlap();
  RealType getMeanOverlap();
  RealType getRelativeVolumeError();

  /** measures over individual labels */
  RealType getSensitivity(LabelType);
  RealType getSpecificity(LabelType);
  RealType getPPV(LabelType);
  RealType getNPV(LabelType);
  RealType getUnionOverlap(LabelType);
  RealType getMeanOverlap(LabelType);
  RealType getRelativeVolumeError(LabelType);

  /** alternative names */
  RealType GetJaccardCoefficient() { return this->getUnionOverlap(); }

  RealType GetDiceCoefficient() { return this->getMeanOverlap(); }

  RealType GetDiceCoefficient(LabelType label) {
    return this->getMeanOverlap(label);
  }

  RealType GetJaccardCoefficient(LabelType label) {
    return this->getUnionOverlap(label);
  }

protected:
  SegmentationMeasuresImageFilter();
  ~SegmentationMeasuresImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void AfterThreadedGenerateData() ITK_OVERRIDE;

  /** Multi-thread version GenerateData. */
  void DynamicThreadedGenerateData(const RegionType &) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(SegmentationMeasuresImageFilter);

  double m_fNbOfPixels;

  std::vector<MapType> m_LabelSetMeasuresPerThread;
  MapType m_LabelSetMeasures;
}; // end of class

} // end namespace anima

#include "animaSegmentationMeasuresImageFilter.hxx"
