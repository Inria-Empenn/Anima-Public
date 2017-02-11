#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkImage.h>

namespace anima
{

template <typename TPixelType, typename TOutputPixelType, unsigned int Dimension>
class JacobianMatrixImageFilter :
public itk::ImageToImageFilter< itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> ,
        itk::Image <itk::Vector <TOutputPixelType, Dimension * Dimension>, Dimension> >
{
public:
    typedef JacobianMatrixImageFilter Self;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> InputImageType;
    typedef typename itk::Image <itk::Vector <TOutputPixelType, Dimension * Dimension>, Dimension> OutputImageType;
    typedef itk::ImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(JacobianMatrixImageFilter, itk::ImageToImageFilter)

    typedef typename InputImageType::IndexType IndexType;
    typedef typename InputImageType::RegionType RegionType;
    typedef typename InputImageType::PointType PointType;
    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction <InputImageType> InterpolatorType;
    typedef typename InterpolatorType::Pointer InterpolatorPointer;

    itkSetMacro(NoIdentity, bool)
    itkSetMacro(Neighborhood, unsigned int)

protected:
    JacobianMatrixImageFilter()
    {
        m_NoIdentity = false;
        m_Neighborhood = 1;
    }

    virtual ~JacobianMatrixImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    JacobianMatrixImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    //! Add identity to the obtained Jacobian matrix
    bool m_NoIdentity;

    //! Neighborhood from which to compute Jacobian matrix
    unsigned int m_Neighborhood;

    InterpolatorPointer m_FieldInterpolator;
};

} // end namespace anima

#include "animaJacobianMatrixImageFilter.hxx"
