#pragma once

namespace anima
{
    /**
     * Constructor
     */
    template <typename TInputImage>
    TensorAsymmetricBlockMatchingRegistrationMethod<TInputImage>
    ::TensorAsymmetricBlockMatchingRegistrationMethod()
    {
    }
    
    template <typename TInputImage>
    void
    TensorAsymmetricBlockMatchingRegistrationMethod<TInputImage>
    ::PerformOneIteration(InputImageType *refImage, InputImageType *movingImage,
                          TransformPointer &addOn)
    {
        // Do a threaded block match computation
        this->GetAgregator()->SetInputRegions(this->GetFixedImage(), this->GetBlockRegions());
        itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
        
        ThreadedMatchData data;
        data.BlockMatch = this;
        data.fixedImage = this->GetFixedImage();
        data.movingImage = movingImage;
        
        itk::TimeProbe tmpTime;
        tmpTime.Start();
        
        unsigned int numThreads = std::min(this->GetNumberOfThreads(),(unsigned int)this->GetBlockRegions().size());

        threader->SetNumberOfThreads(numThreads);
        threader->SetSingleMethod(this->ThreadedMatching, &data);
        threader->SingleMethodExecute();
        
        tmpTime.Stop();
        std::cout << "Matching performed in " << tmpTime.GetTotal() << std::endl;

        this->GetAgregator()->SetInputWeights(this->GetBlockWeights());
        this->GetAgregator()->SetInputTransforms(this->GetBaseTransformsPointers());
        
        addOn = this->GetAgregator()->GetOutput();
    }
    
} // end namespace anima
