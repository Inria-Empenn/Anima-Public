#pragma once

#include <itkImageRegionSplitterBase.h>
#include <vector>

namespace anima
{

template <typename TMaskImage>
class ImageRegionSplitterMask : public itk::ImageRegionSplitterBase
{
public:
    typedef ImageRegionSplitterMask Self;
    typedef itk::ImageRegionSplitterBase Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    typedef TMaskImage MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename TMaskImage::RegionType ImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ImageRegionSplitterMask, ImageRegionSplitterBase);

    void InitializeFromMask(MaskImageType *mask, const ImageRegionType &requestedRegion,
                            unsigned int numberOfThreads);

    unsigned int GetNumberOfThreads()
    {
        return m_lowerLimits.size();
    }

protected:
    ImageRegionSplitterMask()
    {
        m_lowerLimits.clear();
        m_upperLimits.clear();

        m_HighestDim = TMaskImage::ImageDimension;
    }

    virtual unsigned int GetNumberOfSplitsInternal(unsigned int dim,
                                                   const itk::IndexValueType regionIndex[],
                                                   const itk::SizeValueType regionSize[],
                                                   unsigned int requestedNumber ) const;

    virtual unsigned int GetSplitInternal(unsigned int dim,
                                          unsigned int i,
                                          unsigned int numberOfPieces,
                                          itk::IndexValueType regionIndex[],
                                          itk::SizeValueType regionSize[] ) const;

private:
    ImageRegionSplitterMask(const ImageRegionSplitterMask &); //purposely not implemented
    void operator=(const ImageRegionSplitterMask &);      //purposely not implemented

    std::vector <unsigned int> m_lowerLimits, m_upperLimits;

    // Dimension on which the split is done
    unsigned int m_HighestDim;
};

} // end namespace anima

#include "animaImageRegionSplitterMask.hxx"
