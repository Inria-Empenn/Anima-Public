#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

/**
 * @brief Computes the exponentiation of a stationary velocity field using sclaing and squaring
 * and approximated exponential integrators
 *
 * Depending on the order set (0 or 1), the approximated exponentiation for a small enough field
 * is performed using Ferraris et al. ss_aei or Arsigny et al. original approach
 *
 * S. Ferraris et al. Accurate small deformation exponential approximant to integrate large velocity fields: Application to image registration. WBIR 2016
 * V. Arsigny et al. A Log-Euclidean Framework for Statistics on Diffeomorphisms. MICCAI 2006.
 */
template <typename TPixelType, unsigned int Dimension>
class SVFExponentialImageFilter :
public itk::ImageToImageFilter< itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> ,
        itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> >
{
public:
    typedef SVFExponentialImageFilter Self;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> InputImageType;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension * Dimension>, Dimension> JacobianImageType;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> OutputImageType;
    typedef itk::ImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(SVFExponentialImageFilter, itk::ImageToImageFilter)

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;
    typedef typename JacobianImageType::Pointer JacobianImagePointer;
    typedef typename JacobianImageType::PixelType JacobianPixelType;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(ExponentiationOrder, unsigned int)
    itkSetMacro(MaximalDisplacementAmplitude, double)

protected:
    SVFExponentialImageFilter()
    {
        m_ExponentiationOrder = 0;
        m_MaximalDisplacementAmplitude = 0.25;
        m_FieldJacobian = 0;
    }

    virtual ~SVFExponentialImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;
    void AfterThreadedGenerateData() ITK_OVERRIDE;

private:
    SVFExponentialImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    //! Maximal voxel displacement amplitude to be achieved after scaling
    double m_MaximalDisplacementAmplitude;

    //! Exponentiation order (0: Arsigny et al., 1: ss_aei from Ferraris et al.)
    double m_ExponentiationOrder;

    //! Jacobian field (computed only if order 1)
    JacobianImagePointer m_FieldJacobian;

    //! Internal variable that holds the automatically computed number of recursive squarings
    unsigned int m_NumberOfSquarings;
};

} // end namespace anima

#include "animaSVFExponentialImageFilter.hxx"
