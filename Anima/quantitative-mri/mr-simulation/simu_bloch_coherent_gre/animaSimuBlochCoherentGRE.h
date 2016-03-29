#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

template< class TImage>
class SimuBlochCoherentGRE:public itk::ImageToImageFilter< TImage, TImage >
{
public:
    /** Standard class typedefs. */
    typedef SimuBlochCoherentGRE             Self;
    typedef itk::ImageToImageFilter< TImage, TImage > Superclass;
    typedef itk::SmartPointer< Self >        Pointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SimuBlochCoherentGRE, ImageToImageFilter);

    itkSetMacro(TR, float);
    itkGetMacro(TR, float);

    itkSetMacro(TE, float);
    itkGetMacro(TE, float);

    itkSetMacro(FA, float);//changed for CoherentGRE
    itkGetMacro(FA, float);//changed for CoherentGRE


    /** T1 map */
    void SetInputT1(const TImage* T1);

    /** T2s map */
    void SetInputT2s(const TImage* T2s);

    /** M0 image / Rho map */
    void SetInputM0(const TImage* M0);

    /** T2 map */
    void SetInputT2(const TImage* T2);//changed for CoherentGRE

protected:
    SimuBlochCoherentGRE();
    virtual ~SimuBlochCoherentGRE(){}

    /** Does the real work. */
    virtual void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                                      itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    SimuBlochCoherentGRE(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented

    float m_TR;
    float m_TE;
    float m_FA;//changed for CoherentGRE

};

} // end of namespace anima

#include "animaSimuBlochCoherentGRE.hxx"

