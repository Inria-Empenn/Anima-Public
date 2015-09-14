#pragma once
#include "animaRefactoredAsymmetricBMRegistrationMethod.h"

namespace anima
{

template <typename TInputImageType>
void
RefactoredAsymmetricBMRegistrationMethod <TInputImageType>
::PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn)
{
    this->GetAgregator()->SetInputRegions(this->GetFixedImage(), this->GetBlockMatcher()->GetBlockRegions());

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    this->GetBlockMatcher()->SetReferenceImage(refImage);
    this->GetBlockMatcher()->SetMovingImage(movingImage);
    this->GetBlockMatcher()->Update();

    tmpTime.Stop();
    std::cout << "Matching performed in " << tmpTime.GetTotal() << std::endl;

    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());

    addOn = this->GetAgregator()->GetOutput();
}

}
