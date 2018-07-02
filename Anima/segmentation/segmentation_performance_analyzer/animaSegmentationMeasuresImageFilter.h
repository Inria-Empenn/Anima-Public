#pragma once

#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>

namespace anima
{

/** \class LabelSetMeasures
   * \brief Metrics stored per label
   * \ingroup ITKImageStatistics
   */
class SegPerfLabelSetMeasures
{
public:
    // default constructor
    SegPerfLabelSetMeasures()
    {
        m_Source = 0;
        m_Target = 0;
        m_Union = 0;
        m_TrueNegative = 0;
        m_Intersection = 0;
        m_SourceComplement = 0;
        m_TargetComplement = 0;
    }

    // added for completeness
    SegPerfLabelSetMeasures& operator=(const SegPerfLabelSetMeasures& l)
    {
        if(this != &l)
        {
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

template<typename TLabelImage>
class SegmentationMeasuresImageFilter : public itk::ImageToImageFilter<TLabelImage, TLabelImage>
{
public:
    /** Standard Self typedef */
    typedef SegmentationMeasuresImageFilter Self;
    typedef itk::ImageToImageFilter<TLabelImage, TLabelImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(SegmentationMeasuresImageFilter, itk::ImageToImageFilter)

    /** Image related typedefs. */
    typedef TLabelImage LabelImageType;
    typedef typename TLabelImage::Pointer LabelImagePointer;
    typedef typename TLabelImage::ConstPointer LabelImageConstPointer;

    typedef typename TLabelImage::RegionType RegionType;
    typedef typename TLabelImage::SizeType SizeType;
    typedef typename TLabelImage::IndexType IndexType;

    typedef typename TLabelImage::PixelType LabelType;

    /** Type to use for computations. */
    typedef typename itk::NumericTraits<LabelType>::RealType RealType;

    /** Type of the map used to store data per label */
    typedef std::map<LabelType, SegPerfLabelSetMeasures> MapType;
    typedef typename MapType::iterator MapIterator;
    typedef typename MapType::const_iterator MapConstIterator;

    /** Image related typedefs. */
    itkStaticConstMacro(ImageDimension, unsigned int, TLabelImage::ImageDimension);

    /** Set the source image. */
    void SetSourceImage(const LabelImageType* image)
    {
        this->SetNthInput(0, const_cast<LabelImageType*>(image));
    }

    /** Set the target image. */
    void SetTargetImage(const LabelImageType* image)
    {
        this->SetNthInput(1, const_cast<LabelImageType*>(image));
    }

    /** Get the source image. */
    const LabelImageType* GetSourceImage( void )
    {
        return this->GetInput(0);
    }

    /** Get the target image. */
    const LabelImageType* GetTargetImage( void )
    {
        return this->GetInput(1);
    }

    /** Get the label set measures */
    MapType GetLabelSetMeasures()
    {
        return this->m_LabelSetMeasures;
    }

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
    RealType getUnionOverlap( LabelType );
    RealType getMeanOverlap( LabelType );
    RealType getRelativeVolumeError( LabelType );

    /** alternative names */
    RealType GetJaccardCoefficient()
    {
        return this->getUnionOverlap();
    }

    RealType GetDiceCoefficient()
    {
        return this->getMeanOverlap();
    }

    RealType GetDiceCoefficient( LabelType label )
    {
        return this->getMeanOverlap( label );
    }

    RealType GetJaccardCoefficient( LabelType label )
    {
        return this->getUnionOverlap( label );
    }

protected:
    SegmentationMeasuresImageFilter();
    ~SegmentationMeasuresImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void AfterThreadedGenerateData() ITK_OVERRIDE;

    /** Multi-thread version GenerateData. */
    void ThreadedGenerateData(const RegionType &, itk::ThreadIdType) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(SegmentationMeasuresImageFilter);

    float m_fNbOfPixels;

    std::vector<MapType> m_LabelSetMeasuresPerThread;
    MapType m_LabelSetMeasures;
}; // end of class

} // end namespace anima

#include "animaSegmentationMeasuresImageFilter.hxx"
