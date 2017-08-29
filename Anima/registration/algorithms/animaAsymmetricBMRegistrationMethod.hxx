#pragma once
#include "animaAsymmetricBMRegistrationMethod.h"

#include <animaBalooSVFTransformAgregator.h>
#include <animaDenseSVFTransformAgregator.h>

namespace anima
{

template <typename TInputImageType>
void
AsymmetricBMRegistrationMethod <TInputImageType>
::PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn)
{
    itk::TimeProbe tmpTime;
    tmpTime.Start();

    this->GetBlockMatcher()->SetForceComputeBlocks(false);
    this->GetBlockMatcher()->SetReferenceImage(this->GetFixedImage());
    this->GetBlockMatcher()->SetMovingImage(movingImage);
    this->GetBlockMatcher()->SetNumberOfThreads(this->GetNumberOfThreads());
    this->GetBlockMatcher()->Update();

    tmpTime.Stop();

    if (this->GetVerboseProgression())
        std::cout << "Matching performed in " << tmpTime.GetTotal() << std::endl;

    this->GetAgregator()->SetInputRegions(this->GetBlockMatcher()->GetBlockRegions());
    this->GetAgregator()->SetInputOrigins(this->GetBlockMatcher()->GetBlockPositions());

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        typedef anima::BalooSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
        SVFAgregatorType *tmpAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

        if (tmpAgreg)
            tmpAgreg->SetBlockDamWeights(this->GetBlockMatcher()->GetBlockDamWeights());
        else
        {
            typedef anima::DenseSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
            SVFAgregatorType *tmpDenseAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

            tmpDenseAgreg->SetBlockDamWeights(this->GetBlockMatcher()->GetBlockDamWeights());
        }
    }

    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());
    
    addOn = this->GetAgregator()->GetOutput();
}

}
