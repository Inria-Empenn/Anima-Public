#pragma once

#include <itkImageToImageFilter.h>
#include <complex>
#include <vector>

namespace anima
{
    
template <class TImage, class TOutputImage>
class StimulatedSpinEchoImageFilter : public itk::ImageToImageFilter <TImage, TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef StimulatedSpinEchoImageFilter Self;
    typedef itk::ImageToImageFilter <TImage, TOutputImage> Superclass;
    typedef itk::SmartPointer <Self> Pointer;

    typedef itk::Image <typename TImage::PixelType, 4> Image4DType;
    typedef typename Image4DType::Pointer Image4DPointer;
    typedef TOutputImage OutputImageType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef std::vector < std::complex <double> > ComplexVectorType;
    typedef std::vector <ComplexVectorType> MatrixType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(StimulatedSpinEchoImageFilter, itk::ImageToImageFilter)

    itkSetMacro(EchoSpacing, double)
    itkGetMacro(EchoSpacing, double)

    itkSetMacro(NumberOfEchoes, unsigned int)
    itkGetMacro(NumberOfEchoes, unsigned int)

    itkSetMacro(FlipAngle, double)
    itkGetMacro(FlipAngle, double)

    itkSetMacro(ExcitationFlipAngle, double)
    itkGetMacro(ExcitationFlipAngle, double)

    itkSetMacro(B1OnExcitationAngle, bool)

    /** T1 map */
    void SetInputT1(const TImage* T1);

    /** T2 map */
    void SetInputT2(const TImage* T2);

    /** M0 image / Rho map */
    void SetInputM0(const TImage* M0);

    /** B1 inhomogeneity image */
    void SetInputB1(const TImage* B1);

    Image4DType *GetOutputAs4DImage();

protected:
    StimulatedSpinEchoImageFilter();
    virtual ~StimulatedSpinEchoImageFilter() {}

    /** Does the real work. */
    virtual void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                                      itk::ThreadIdType threadId) ITK_OVERRIDE;

    void GenerateOutputInformation() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(StimulatedSpinEchoImageFilter);

    double m_EchoSpacing;
    double m_ExcitationFlipAngle;
    double m_FlipAngle;
    unsigned int m_NumberOfEchoes;

    bool m_B1OnExcitationAngle;

    Image4DPointer m_Output4D;
};
    
} // end of namespace anima

#include "animaStimulatedSpinEchoImageFilter.hxx"
