#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

/**
 * @brief Computes the Lie bracket between two fields u and v as expressed by Bossa et al.
 *
 * There is a discrepancy in between Vercauteren et al. and Bossa et al. formulation. This Lie bracket implements the
 * Vercauteren et al. formulation
 * [u,v](x) = Jac(u)(x).v(x) - Jac(v)(x).u(x)
 * M. Bossa et al. "Contributions to 3D diffeomorphic atlas estimation : application to brain images.", MICCAI 2007, p. 667–674.
 * T. Vercauteren et al. "Symmetric Log-Domain Diffeomorphic Registration: A Demons-based Approach.", MICCAI 2008, p. 754-761.
 */
template <typename TPixelType, unsigned int Dimension>
class SVFLieBracketImageFilter :
public itk::ImageToImageFilter< itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> ,
        itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> >
{
public:
    typedef SVFLieBracketImageFilter Self;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> InputImageType;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension * Dimension>, Dimension> JacobianImageType;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> OutputImageType;
    typedef itk::ImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(SVFLieBracketImageFilter, itk::ImageToImageFilter)

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;
    typedef typename JacobianImageType::Pointer JacobianImagePointer;
    typedef typename JacobianImageType::PixelType JacobianPixelType;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetObjectMacro(FirstFieldJacobian, JacobianImageType)
    itkSetObjectMacro(SecondFieldJacobian, JacobianImageType)
    itkGetObjectMacro(FirstFieldJacobian, JacobianImageType)
    itkGetObjectMacro(SecondFieldJacobian, JacobianImageType)

protected:
    SVFLieBracketImageFilter()
    {
        this->SetNumberOfRequiredInputs(2);
        m_FirstFieldJacobian = 0;
        m_SecondFieldJacobian = 0;
    }

    virtual ~SVFLieBracketImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    SVFLieBracketImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    JacobianImagePointer m_FirstFieldJacobian, m_SecondFieldJacobian;
};

} // end namespace anima

#include "animaSVFLieBracketImageFilter.hxx"
