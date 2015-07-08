#pragma once

#include <itkTimeProbe.h>

#include <animaMatrixLogExp.h>

#include <animaBalooSVFTransformAgregator.h>
#include <animaDenseSVFTransformAgregator.h>

namespace anima
{
/**
 * Constructor
 */
template <typename TInputImage>
TensorKissingSymmetricBlockMatchingRegistrationMethod<TInputImage>
::TensorKissingSymmetricBlockMatchingRegistrationMethod()
{
    m_BlockPercentageKept = 0.8;
    m_BlockSize = 5;
    m_BlockSpacing = 2;
    m_UseTransformationDam = true;
    m_DamDistance = 6;
}

template <typename TInputImage>
void
TensorKissingSymmetricBlockMatchingRegistrationMethod<TInputImage>
::InitializeBlocksOnImage(InitializerPointer &initPtr, InputImageType *image)
{
    // Init blocks
    initPtr = InitializerType::New();
    initPtr->AddReferenceImage(image);

    if (this->GetNumberOfThreads() != 0)
        initPtr->SetNumberOfThreads(this->GetNumberOfThreads());

    initPtr->SetPercentageKept(m_BlockPercentageKept);
    initPtr->SetBlockSize(m_BlockSize);
    initPtr->SetBlockSpacing(m_BlockSpacing);
    initPtr->SetTensorVarianceThreshold(this->GetBlockTensorVarianceThreshold());

    initPtr->SetRequestedRegion(image->GetLargestPossibleRegion());

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        initPtr->SetComputeOuterDam(m_UseTransformationDam);
        initPtr->SetDamDistance(m_DamDistance);
    }
}

template <typename TInputImage>
void
TensorKissingSymmetricBlockMatchingRegistrationMethod<TInputImage>
::PerformOneIteration(InputImageType *refImage, InputImageType *movingImage,
                      TransformPointer &addOn)
{
    // First initialize blocks (and dams if needed)
    InitializerPointer initPtr;
    this->InitializeBlocksOnImage(initPtr,refImage);

    std::vector <ImageRegionType> fixedBlockRegions = initPtr->GetOutput();
    std::vector <ImageIndexType> fixedDamIndexes = initPtr->GetDamIndexes();

    std::cout << "Generated " << fixedBlockRegions.size() << " blocks..." << std::endl;

    // Do a threaded block match computation for forward matching
    this->SetBlockRegions(fixedBlockRegions);
    this->GetAgregator()->SetInputRegions(refImage, fixedBlockRegions);

    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

    ThreadedMatchData data;
    data.BlockMatch = this;
    data.fixedImage = refImage;
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

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        typedef BalooSVFTransformAgregator<TInputImage::ImageDimension> SVFAgregatorType;
        SVFAgregatorType *tmpAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

        if (tmpAgreg)
            tmpAgreg->SetDamIndexes(fixedDamIndexes);
        else
        {
            typedef DenseSVFTransformAgregator<TInputImage::ImageDimension> SVFAgregatorType;
            SVFAgregatorType *tmpDenseAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

            tmpDenseAgreg->SetDamIndexes(fixedDamIndexes);
        }
    }

    TransformPointer usualAddOn = this->GetAgregator()->GetOutput();

    // Initialize blocks (and dams if needed)
    this->InitializeBlocksOnImage(initPtr,movingImage);

    std::vector <ImageRegionType> movingBlockRegions = initPtr->GetOutput();
    std::vector <ImageIndexType> movingDamIndexes = initPtr->GetDamIndexes();

    std::cout << "Generated " << movingBlockRegions.size() << " blocks..." << std::endl;

    // Do a threaded block match computation for backward matching
    this->SetBlockRegions(movingBlockRegions);
    this->GetAgregator()->SetInputRegions(movingImage, movingBlockRegions);

    threader = itk::MultiThreader::New();

    data.BlockMatch = this;
    data.fixedImage = movingImage;
    data.movingImage = refImage;

    itk::TimeProbe tmpTimeReverse;
    tmpTimeReverse.Start();

    threader->SetNumberOfThreads(numThreads);
    threader->SetSingleMethod(this->ThreadedMatching, &data);
    threader->SingleMethodExecute();

    tmpTimeReverse.Stop();
    std::cout << "Matching performed in " << tmpTimeReverse.GetTotal() << std::endl;

    this->GetAgregator()->SetInputWeights(this->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBaseTransformsPointers());

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        typedef BalooSVFTransformAgregator<TInputImage::ImageDimension> SVFAgregatorType;
        SVFAgregatorType *tmpAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

        if (tmpAgreg)
            tmpAgreg->SetDamIndexes(movingDamIndexes);
        else
        {
            typedef DenseSVFTransformAgregator<TInputImage::ImageDimension> SVFAgregatorType;
            SVFAgregatorType *tmpDenseAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

            tmpDenseAgreg->SetDamIndexes(movingDamIndexes);
        }
    }

    TransformPointer reverseAddOn = this->GetAgregator()->GetOutput();

    if (this->GetAgregator()->GetOutputTransformType() == AgregatorType::SVF)
    {
        // Add update to current velocity field (cf. Vercauteren et al, 2008)
        // First compute the SVF from two asymmetric ones: S = 0.25 * (S_0 - S_1)
        // It's only a quarter since we are computing the half power of the transform between the two images
        typedef typename SVFTransformType::VectorFieldType VelocityFieldType;

        typedef typename itk::ImageRegionConstIterator <VelocityFieldType> VelocityFieldConstIterator;
        typedef typename itk::ImageRegionIterator <VelocityFieldType> VelocityFieldIterator;

        SVFTransformType *usualAddOnCast = dynamic_cast <SVFTransformType *> (usualAddOn.GetPointer());
        SVFTransformType *reverseAddOnCast = dynamic_cast <SVFTransformType *> (reverseAddOn.GetPointer());

        VelocityFieldConstIterator usualAddOnItr(usualAddOnCast->GetParametersAsVectorField(),
                                                 usualAddOnCast->GetParametersAsVectorField()->GetLargestPossibleRegion());

        VelocityFieldConstIterator reverseAddOnItr(reverseAddOnCast->GetParametersAsVectorField(),
                                                   reverseAddOnCast->GetParametersAsVectorField()->GetLargestPossibleRegion());

        typename VelocityFieldType::Pointer resField = VelocityFieldType::New();
        resField->Initialize();
        resField->SetRegions (usualAddOnCast->GetParametersAsVectorField()->GetLargestPossibleRegion());
        resField->SetSpacing (usualAddOnCast->GetParametersAsVectorField()->GetSpacing());
        resField->SetOrigin (usualAddOnCast->GetParametersAsVectorField()->GetOrigin());
        resField->SetDirection (usualAddOnCast->GetParametersAsVectorField()->GetDirection());
        resField->Allocate();

        VelocityFieldIterator resItr(resField,resField->GetLargestPossibleRegion());

        typedef typename VelocityFieldType::PixelType VectorType;
        VectorType tmpVec;
        while (!usualAddOnItr.IsAtEnd())
        {
            tmpVec = 0.25 * (usualAddOnItr.Get() - reverseAddOnItr.Get());
            resItr.Set(tmpVec);

            ++usualAddOnItr;
            ++reverseAddOnItr;
            ++resItr;
        }

        usualAddOnCast->SetParametersAsVectorField(resField.GetPointer());
    }
    else
    {
        AffineTransformType *usualAddOnCast = dynamic_cast <AffineTransformType *> (usualAddOn.GetPointer());
        AffineTransformType *reverseAddOnCast = dynamic_cast <AffineTransformType *> (reverseAddOn.GetPointer());

        unsigned int NDimensions = TInputImage::GetImageDimension();
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

} // end of namespace anima
