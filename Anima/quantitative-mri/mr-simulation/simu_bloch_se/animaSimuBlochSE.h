#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

template< class TImage>
class SimuBlochSE:public itk::ImageToImageFilter< TImage, TImage >
{
public:
    /** Standard class typedefs. */
    typedef SimuBlochSE Self;
    typedef itk::ImageToImageFilter <TImage, TImage> Superclass;
    typedef itk::SmartPointer <Self> Pointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(SimuBlochSE, ImageToImageFilter)

    itkSetMacro(TR, float)
    itkGetMacro(TR, float)

    itkSetMacro(TE, float)
    itkGetMacro(TE, float)

    /** T1 map */
    void SetInputT1(const TImage* T1);

    /** T2 map */
    void SetInputT2(const TImage* T2);

    /** M0 image / Rho map */
    void SetInputM0(const TImage* M0);

protected:
    SimuBlochSE();
    virtual ~SimuBlochSE() {}

    /** Does the real work. */
    virtual void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(SimuBlochSE);

    float m_TR;
    float m_TE;
};

} // end of namespace anima

#include "animaSimuBlochSE.hxx"
