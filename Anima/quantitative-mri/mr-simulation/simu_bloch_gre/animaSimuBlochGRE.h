#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

template< class TImage>
class SimuBlochGRE:public itk::ImageToImageFilter< TImage, TImage >
{
public:
    /** Standard class typedefs. */
    typedef SimuBlochGRE Self;
    typedef itk::ImageToImageFilter <TImage, TImage> Superclass;
    typedef itk::SmartPointer <Self> Pointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(SimuBlochGRE, ImageToImageFilter)

    itkSetMacro(TR, double)
    itkGetMacro(TR, double)

    itkSetMacro(TE, double)
    itkGetMacro(TE, double)

    /** T1 map */
    void SetInputT1(const TImage* T1);

    /** T2s map */
    void SetInputT2s(const TImage* T2s);//changed for GRE

    /** M0 image / Rho map */
    void SetInputM0(const TImage* M0);

protected:
    SimuBlochGRE();
    virtual ~SimuBlochGRE(){}

    /** Does the real work. */
    virtual void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(SimuBlochGRE);

    double m_TR;
    double m_TE;
};

} // end of namespace anima

#include "animaSimuBlochGRE.hxx"
