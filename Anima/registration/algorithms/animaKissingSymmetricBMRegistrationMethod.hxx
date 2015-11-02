#pragma once
#include "animaKissingSymmetricBMRegistrationMethod.h"

#include <animaBalooSVFTransformAgregator.h>
#include <animaDenseSVFTransformAgregator.h>

#include <itkMultiplyImageFilter.h>
#include <itkSubtractImageFilter.h>

namespace anima
{

template <typename TInputImageType>
void
KissingSymmetricBMRegistrationMethod <TInputImageType>
::PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn)
{
    itk::TimeProbe tmpTime;
    tmpTime.Start();

    this->GetBlockMatcher()->SetForceComputeBlocks(true);
    this->GetBlockMatcher()->SetReferenceImage(refImage);
    this->GetBlockMatcher()->SetMovingImage(movingImage);
    this->GetBlockMatcher()->SetNumberOfThreads(this->GetNumberOfThreads());
    this->GetBlockMatcher()->Update();

    tmpTime.Stop();
    std::cout << "Matching performed in " << tmpTime.GetTotal() << std::endl;

    this->GetAgregator()->SetInputRegions(this->GetBlockMatcher()->GetBlockRegions());
    this->GetAgregator()->SetInputOrigins(this->GetBlockMatcher()->GetBlockPositions());

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        typedef anima::BalooSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
        SVFAgregatorType *tmpAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

        if (tmpAgreg)
            tmpAgreg->SetDamIndexes(this->GetBlockMatcher()->GetDamIndexes());
        else
        {
            typedef anima::DenseSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
            SVFAgregatorType *tmpDenseAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

            tmpDenseAgreg->SetDamIndexes(this->GetBlockMatcher()->GetDamIndexes());
        }
    }

    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());

    TransformPointer usualAddOn = this->GetAgregator()->GetOutput();

    itk::TimeProbe tmpTimeReverse;
    tmpTimeReverse.Start();

    this->GetBlockMatcher()->SetReferenceImage(movingImage);
    this->GetBlockMatcher()->SetMovingImage(refImage);
    this->GetBlockMatcher()->Update();

    tmpTimeReverse.Stop();
    std::cout << "Matching performed in " << tmpTimeReverse.GetTotal() << std::endl;

    this->GetAgregator()->SetInputRegions(this->GetBlockMatcher()->GetBlockRegions());
    this->GetAgregator()->SetInputOrigins(this->GetBlockMatcher()->GetBlockPositions());

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        typedef anima::BalooSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
        SVFAgregatorType *tmpAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

        if (tmpAgreg)
            tmpAgreg->SetDamIndexes(this->GetBlockMatcher()->GetDamIndexes());
        else
        {
            typedef anima::DenseSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
            SVFAgregatorType *tmpDenseAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

            tmpDenseAgreg->SetDamIndexes(this->GetBlockMatcher()->GetDamIndexes());
        }
    }

    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());

    TransformPointer reverseAddOn = this->GetAgregator()->GetOutput();

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        // Add update to current velocity field (cf. Vercauteren et al, 2008)
        // First compute the SVF from two asymmetric ones: S = 0.25 * (S_0 - S_1)
        // It's only a quarter since we are computing the half power of the transform between the two images
        typedef typename SVFTransformType::VectorFieldType VelocityFieldType;

        typedef itk::MultiplyImageFilter <VelocityFieldType,itk::Image <float,InputImageType::ImageDimension>,VelocityFieldType> MultiplyFilterType;
        typedef itk::SubtractImageFilter <VelocityFieldType,VelocityFieldType,VelocityFieldType> SubtractFilterType;

        SVFTransformType *usualAddOnCast = dynamic_cast <SVFTransformType *> (usualAddOn.GetPointer());
        SVFTransformType *reverseAddOnCast = dynamic_cast <SVFTransformType *> (reverseAddOn.GetPointer());

        typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
        subFilter->SetInput1(usualAddOnCast->GetParametersAsVectorField());
        subFilter->SetInput2(reverseAddOnCast->GetParametersAsVectorField());
        subFilter->SetNumberOfThreads(this->GetNumberOfThreads());
        subFilter->InPlaceOn();

        subFilter->Update();

        typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
        multiplyFilter->SetInput(subFilter->GetOutput());
        multiplyFilter->SetConstant(0.25);
        multiplyFilter->SetNumberOfThreads(this->GetNumberOfThreads());
        multiplyFilter->InPlaceOn();

        multiplyFilter->Update();

        VelocityFieldType *addonSVF = multiplyFilter->GetOutput();
        addonSVF->DisconnectPipeline();

        usualAddOnCast->SetParametersAsVectorField(addonSVF);
    }
    else
    {
        AffineTransformType *usualAddOnCast = dynamic_cast <AffineTransformType *> (usualAddOn.GetPointer());
        AffineTransformType *reverseAddOnCast = dynamic_cast <AffineTransformType *> (reverseAddOn.GetPointer());

        unsigned int NDimensions = InputImageType::ImageDimension;
        vnl_matrix <double> usualAddOnMatrix(NDimensions+1,NDimensions+1,0);
        vnl_matrix <double> reverseAddOnMatrix(NDimensions+1,NDimensions+1,0);
        usualAddOnMatrix.set_identity();
        reverseAddOnMatrix.set_identity();

        for (unsigned int i = 0;i < NDimensions;++i)
        {
            for (unsigned int j = 0;j < NDimensions;++j)
            {
                usualAddOnMatrix(i,j) = usualAddOnCast->GetMatrix()(i,j);
                reverseAddOnMatrix(i,j) = reverseAddOnCast->GetMatrix()(i,j);
            }

            usualAddOnMatrix(i,NDimensions) = usualAddOnCast->GetOffset()[i];
            reverseAddOnMatrix(i,NDimensions) = reverseAddOnCast->GetOffset()[i];
        }

        usualAddOnMatrix = anima::GetLogarithm(usualAddOnMatrix);
        reverseAddOnMatrix = anima::GetLogarithm(reverseAddOnMatrix);

        usualAddOnMatrix -= reverseAddOnMatrix;
        usualAddOnMatrix /= 4.0;

        usualAddOnMatrix = anima::GetExponential(usualAddOnMatrix);

        typename AffineTransformType::MatrixType trsfMatrix;

        for (unsigned int i = 0;i < NDimensions;++i)
            for (unsigned int j = 0;j < NDimensions;++j)
                trsfMatrix(i,j) = usualAddOnMatrix(i,j);

        usualAddOnCast->SetMatrix(trsfMatrix);

        typename AffineTransformType::OffsetType trsfOffset;
        for (unsigned int i = 0;i < NDimensions;++i)
            trsfOffset[i] = usualAddOnMatrix(i,NDimensions);

        usualAddOnCast->SetOffset(trsfOffset);
    }

    addOn = usualAddOn;
}

}
