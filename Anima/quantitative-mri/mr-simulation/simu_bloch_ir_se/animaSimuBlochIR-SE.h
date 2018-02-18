#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

template< class TImage>
class SimuBlochIRSE:public itk::ImageToImageFilter< TImage, TImage >
{
public:
    /** Standard class typedefs. */
    typedef SimuBlochIRSE Self;
    typedef itk::ImageToImageFilter <TImage, TImage> Superclass;
    typedef itk::SmartPointer <Self> Pointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(SimuBlochIRSE, ImageToImageFilter)

    itkSetMacro(TR, float)
    itkGetMacro(TR, float)

    itkSetMacro(TE, float)
    itkGetMacro(TE, float)

    itkSetMacro(TI, float)
    itkGetMacro(TI, float)

    /** T1 map */
    void SetInputT1(const TImage* T1);

    /** T2 map */
    void SetInputT2(const TImage* T2);

    /** M0 image / Rho map */
    void SetInputM0(const TImage* M0);

protected:
    SimuBlochIRSE();
    virtual ~SimuBlochIRSE(){}

    /** Does the real work. */
    virtual void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                                      itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(SimuBlochIRSE);

    float m_TR;
    float m_TE;
    float m_TI;
};

} // end of namespace anima

#include "animaSimuBlochIR-SE.hxx"
