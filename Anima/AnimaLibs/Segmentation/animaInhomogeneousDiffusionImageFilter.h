#pragma once

#include <animaInhomogeneousAOSLineDiffusionImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkImage.h>
#include <itkPixelTraits.h>
#include <itkCommand.h>

namespace anima
{

template <typename TInputImage, typename TDiffusionScalarImage=TInputImage, typename TOutputImage=TInputImage>
class InhomogeneousDiffusionImageFilter :
        public itk::ImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef InhomogeneousDiffusionImageFilter          Self;
    typedef itk::ImageToImageFilter<TInputImage,TOutputImage>      Superclass;
    typedef itk::SmartPointer<Self>                                Pointer;
    typedef itk::SmartPointer<const Self>                          ConstPointer;

    /** Pixel Type of the input image */
    typedef TInputImage                                       InputImageType;
    typedef TDiffusionScalarImage                             DiffusionScalarsImageType;
    typedef TOutputImage                                      OutputImageType;
    typedef typename TInputImage::PixelType                   PixelType;
    typedef typename itk::NumericTraits<PixelType>::RealType       RealType;
    typedef typename itk::NumericTraits<PixelType>::ScalarRealType ScalarRealType;


    /** Runtime information support. */
    itkTypeMacro(InhomogeneousDiffusionImageFilter,
                 ImageToImageFilter)

    /** Image dimension. */
    itkStaticConstMacro(ImageDimension, unsigned int,
                        TInputImage::ImageDimension);

    /** Define the image type for internal computations
         RealType is usually 'double' in NumericTraits.
         Here we prefer double in order to save memory.  */

    typedef typename itk::NumericTraits< PixelType >::RealType   InternalRealType;
    typedef typename InputImageType::template Rebind<InternalRealType>::Type RealImageType;

    /**  The first in the pipeline  */
    typedef anima::InhomogeneousAOSLineDiffusionImageFilter <
    RealImageType,
    DiffusionScalarsImageType,
    RealImageType
    > InternalAOSDiffusionFilterType;

    typedef itk::CastImageFilter <InputImageType,RealImageType> RealCastingFilterType;
    typedef itk::CastImageFilter <RealImageType,OutputImageType> OutCastingFilterType;

    /**  Pointer to a gaussian filter.  */
    typedef typename InternalAOSDiffusionFilterType::Pointer    InternalAOSDiffusionFilterPointer;

    typedef typename DiffusionScalarsImageType::Pointer     DiffusionScalarsImagePointer;
    typedef typename OutputImageType::Pointer               OutputImagePointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    itkSetMacro(NumberOfSteps, unsigned int)
    itkSetMacro(StepLength, double)
    itkSetMacro(DiffusionSourceFactor, double)
    itkSetObjectMacro(DiffusionScalarsImage, DiffusionScalarsImageType)

protected:
    InhomogeneousDiffusionImageFilter();
    virtual ~InhomogeneousDiffusionImageFilter() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    /** Generate Data */
    void GenerateData() ITK_OVERRIDE;

    /** InhomogeneousDiffusionImageFilter needs all of the input to produce an
         * output. Therefore, InhomogeneousDiffusionImageFilter needs to provide
         * an implementation for GenerateInputRequestedRegion in order to inform
         * the pipeline execution model.
         * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
    virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

    // Override since the filter produces the entire dataset
    void EnlargeOutputRequestedRegion(itk::DataObject *output) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(InhomogeneousDiffusionImageFilter);

    DiffusionScalarsImagePointer m_DiffusionScalarsImage;

    unsigned int m_NumberOfSteps;
    double m_StepLength;
    double m_DiffusionSourceFactor;
};

} // end namespace anima

#include "animaInhomogeneousDiffusionImageFilter.hxx"
