#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

template <class TPixelType>
class MajorityLabelVotingImageFilter :
public itk::ImageToImageFilter< itk::Image <TPixelType, 3>, itk::Image<TPixelType, 3> >
{
public:
    /** Standard class type def */

    using Self = MajorityLabelVotingImageFilter;
    using InputImageType = itk::Image <TPixelType, 3>;
    using OutputImageType = itk::Image <TPixelType, 3>;
    using Superclass = itk::ImageToImageFilter <InputImageType, OutputImageType>;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MajorityLabelVotingImageFilter, ImageToImageFilter)

    using InputRegionType = typename InputImageType::RegionType;

protected:
    MajorityLabelVotingImageFilter ()
    {
    }

    virtual ~MajorityLabelVotingImageFilter () {}

    void DynamicThreadedGenerateData(const InputRegionType &region) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MajorityLabelVotingImageFilter);
};

} // end namespace anima

#include "animaMajorityLabelVotingImageFilter.hxx"
