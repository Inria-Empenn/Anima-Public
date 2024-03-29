#pragma once
#include "animaKissingSymmetricBMRegistrationMethod.h"

#include <animaBalooSVFTransformAgregator.h>
#include <animaDenseSVFTransformAgregator.h>
#include <animaMatrixLogExp.h>

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

    this->ChangeDefaultBackgroundValue(refImage, m_FloatingBackgroundValue);
    this->GetBlockMatcher()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    this->GetBlockMatcher()->Update();

    tmpTime.Stop();

    if (this->GetVerboseProgression())
        std::cout << "Matching performed in " << tmpTime.GetTotal() << std::endl;

    this->GetAgregator()->SetInputRegions(this->GetBlockMatcher()->GetBlockRegions());
    this->GetAgregator()->SetInputOrigins(this->GetBlockMatcher()->GetBlockPositions());
    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());

    TransformPointer usualAddOn = this->GetAgregator()->GetOutput();

    itk::TimeProbe tmpTimeReverse;
    tmpTimeReverse.Start();

    this->GetBlockMatcher()->SetReferenceImage(movingImage);
    this->GetBlockMatcher()->SetMovingImage(refImage);

    this->ChangeDefaultBackgroundValue(refImage, m_ReferenceBackgroundValue);
    this->GetBlockMatcher()->Update();

    tmpTimeReverse.Stop();

    if (this->GetVerboseProgression())
        std::cout << "Matching performed in " << tmpTimeReverse.GetTotal() << std::endl;

    this->GetAgregator()->SetInputRegions(this->GetBlockMatcher()->GetBlockRegions());
    this->GetAgregator()->SetInputOrigins(this->GetBlockMatcher()->GetBlockPositions());
    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());

    TransformPointer reverseAddOn = this->GetAgregator()->GetOutput();

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        // Add update to current velocity field (cf. Vercauteren et al, 2008)
        // First compute the SVF from two asymmetric ones: S = RPL * (S_0 - S_1) / 2
        // It's only a quarter since we are computing the half power of the transform between the two images
        typedef typename SVFTransformType::VectorFieldType VelocityFieldType;

        typedef itk::MultiplyImageFilter <VelocityFieldType,itk::Image <double,InputImageType::ImageDimension>,VelocityFieldType> MultiplyFilterType;
        typedef itk::SubtractImageFilter <VelocityFieldType,VelocityFieldType,VelocityFieldType> SubtractFilterType;

        SVFTransformType *usualAddOnCast = dynamic_cast <SVFTransformType *> (usualAddOn.GetPointer());
        SVFTransformType *reverseAddOnCast = dynamic_cast <SVFTransformType *> (reverseAddOn.GetPointer());

        typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
        subFilter->SetInput1(usualAddOnCast->GetParametersAsVectorField());
        subFilter->SetInput2(reverseAddOnCast->GetParametersAsVectorField());
        subFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        subFilter->InPlaceOn();

        subFilter->Update();

        typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
        multiplyFilter->SetInput(subFilter->GetOutput());
        multiplyFilter->SetConstant(0.5 * m_RegistrationPointLocation);
        multiplyFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        multiplyFilter->InPlaceOn();

        multiplyFilter->Update();

        typename VelocityFieldType::Pointer addonSVF = multiplyFilter->GetOutput();
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
        usualAddOnMatrix *= m_RegistrationPointLocation / 2.0;

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

template <typename TInputImageType>
typename KissingSymmetricBMRegistrationMethod <TInputImageType>::TransformPointer
KissingSymmetricBMRegistrationMethod <TInputImageType>
::GetBackwardTransformForResampling(TransformType *transform)
{
    TransformPointer outTrsf;

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        // Compute temporary field and set it to resampler
        DisplacementFieldTransformPointer dispTrsf = DisplacementFieldTransformType::New();
        SVFTransformType *svfCast = dynamic_cast<SVFTransformType *> (transform);

        anima::GetSVFExponential(svfCast,dispTrsf.GetPointer(),this->GetExponentiationOrder(),this->GetNumberOfWorkUnits(), (m_RegistrationPointLocation - 1) / m_RegistrationPointLocation);

        outTrsf = dispTrsf;
    }
    else
    {
        AffineTransformType *affCast = dynamic_cast<AffineTransformType *> (transform);

        unsigned int NDimensions = InputImageType::ImageDimension;
        vnl_matrix <double> affMatrix(NDimensions+1,NDimensions+1,0);
        affMatrix.set_identity();

        for (unsigned int i = 0;i < NDimensions;++i)
        {
            for (unsigned int j = 0;j < NDimensions;++j)
                affMatrix(i,j) = affCast->GetMatrix()(i,j);

            affMatrix(i,NDimensions) = affCast->GetOffset()[i];
        }

        affMatrix = anima::GetLogarithm(affMatrix);
        affMatrix *= (m_RegistrationPointLocation - 1) / m_RegistrationPointLocation;
        affMatrix = anima::GetExponential(affMatrix);

        vnl_matrix <double> affResult(NDimensions,NDimensions,0);
        typename AffineTransformType::OffsetType offset(NDimensions);

        for (unsigned int i = 0;i < NDimensions;++i)
        {
            for (unsigned int j = 0;j < NDimensions;++j)
                affResult(i,j) = affMatrix(i,j);

            offset[i] = affMatrix(i,NDimensions);
        }

        AffineTransformPointer outCast = AffineTransformType::New();
        outCast->SetMatrix(affResult);
        outCast->SetOffset(offset);

        outTrsf = outCast;
    }

    return outTrsf;
}

}
