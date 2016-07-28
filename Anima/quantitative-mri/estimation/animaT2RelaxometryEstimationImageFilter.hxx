#pragma once
#include "animaT2RelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <nlopt.hpp>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
void
T2RelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;

    typedef itk::ImageRegionConstIterator <InputImageType> IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    std::vector <IteratorType> inItrs(this->GetNumberOfIndexedInputs());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        inItrs[i] = IteratorType(this->GetInput(i),this->GetOutput()->GetLargestPossibleRegion());

    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    maskImage->SetSpacing (this->GetInput(0)->GetSpacing());
    maskImage->SetOrigin (this->GetInput(0)->GetOrigin());
    maskImage->SetDirection (this->GetInput(0)->GetDirection());
    maskImage->Allocate();

    MaskIteratorType maskItr (maskImage,this->GetOutput()->GetLargestPossibleRegion());
    while (!maskItr.IsAtEnd())
    {
        double averageVal = 0;
        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            averageVal += inItrs[i].Get();

        averageVal /= this->GetNumberOfIndexedInputs();

        bool maskPoint = (averageVal <= m_AverageSignalThreshold);

        if (maskPoint)
            maskItr.Set(0);
        else
            maskItr.Set(1);

        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            ++inItrs[i];

        ++maskItr;
    }

    this->SetComputationMask(maskImage);
}

template <typename TInputImage, typename TOutputImage>
void
T2RelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;

    OutImageIteratorType outT2Iterator(this->GetOutput(0),outputRegionForThread);
    OutImageIteratorType outM0Iterator(this->GetOutput(1),outputRegionForThread);

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    OutImageIteratorType t1MapItr;
    if (m_T1Map)
        t1MapItr = OutImageIteratorType(m_T1Map,outputRegionForThread);

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            outT2Iterator.Set(0);
            outM0Iterator.Set(0);

            ++maskItr;
            ++outT2Iterator;
            ++outM0Iterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            continue;
        }

        double t1Value = 1;

        if (m_T1Map)
            t1Value = t1MapItr.Get();


        // must be computed here
        nlopt::opt opt(nlopt::LN_BOBYQA, 2);

        double* echoTimeAndObservation=new double[numInputs*2+1];
        echoTimeAndObservation[0]=numInputs;
        for(unsigned int index=0;index<numInputs;++index)
        {
           echoTimeAndObservation[2*index+1]=this->m_EchoTime[index];
           echoTimeAndObservation[2*index+2]=inIterators[index].Get();
         }

        std::vector<double>lowerBounds(2,1e-4);
        opt.set_lower_bounds(lowerBounds);
        std::vector<double>upperBounds(2);upperBounds[0]=m_T2UpperBoundValue;upperBounds[1]=m_M0UpperBoundValue;
        opt.set_upper_bounds(upperBounds);
        // one must also set upper bounds
        opt.set_min_objective(T2RelaxometryEstimationImageFilter <TInputImage,TOutputImage>::myfunc, echoTimeAndObservation);
        opt.set_xtol_rel(1e-4);

        std::vector<double> optimizedValue(2);
        optimizedValue[0] = 100; optimizedValue[1] = 600;
        double minf;

        double t2Value=optimizedValue[0];
        double m0Value=optimizedValue[1];
        try
        {
            opt.optimize(optimizedValue, minf);
        }
        catch(nlopt::roundoff_limited& e)
        {
            std::cerr << "NLOPT optimization error, at: "<< inIterators[0].GetIndex() << std::endl;
            for (unsigned int i = 0;i < numInputs;++i)
                std::cerr << "inIterators[" <<  i << "]: " << inIterators[i].Get() << std::endl;

            outT2Iterator.Set(lowerBounds[0]);
            outM0Iterator.Set(lowerBounds[1]);

            ++maskItr;
            ++outT2Iterator;
            ++outM0Iterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            continue;
        }

        t2Value=optimizedValue[0];
        m0Value=optimizedValue[1]/(1-exp(-m_TRValue/t1Value));

        outM0Iterator.Set(m0Value);
        outT2Iterator.Set(t2Value);

        ++maskItr;
        ++outT2Iterator;
        ++outM0Iterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;
    }
}

template <typename TInputImage, typename TOutputImage>
double
T2RelaxometryEstimationImageFilter <TInputImage,TOutputImage>::myfunc(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    double * my_func_data=static_cast<double*>(data);
    grad.resize(0);
    if (grad.size()==2) {
        grad[0] = 0.0;
        grad[1] = 0.5 / std::sqrt(x[1]);
    }
    unsigned int nbOfEchos=my_func_data[0];//
    double res=0;
    const double T2=x[0];
    const double scale=x[1];
    for(unsigned int index=0;index<nbOfEchos;++index)
    {
        const double echo=my_func_data[2*index+1];
        const double simulatedSignal=scale*exp(-echo/T2);
        const double observedSignal=my_func_data[2*index+2];
        //std::cout<<observedSignal<<" "<<simulatedSignal<<std::endl;
        res+=(simulatedSignal-observedSignal)*(simulatedSignal-observedSignal);
    }
    //std::cout<<res<<std::endl;
    return res;
}

} // end namespace anima
