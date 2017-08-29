#pragma once
#include "animaBaseBMRegistrationMethod.h"

#include <animaResampleImageFilter.h>
#include <animaOrientedModelBaseResampleImageFilter.h>
#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>

#include <animaVelocityUtils.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <typename TInputImageType>
BaseBMRegistrationMethod <TInputImageType>
::BaseBMRegistrationMethod()
{
    m_Abort = false;
    m_FixedImage = 0;
    m_MovingImage = 0;

    m_SVFElasticRegSigma = 0;
    m_BCHCompositionOrder = 1;
    m_ExponentiationOrder = 0;

    m_MaximumIterations = 10;
    m_MinimalTransformError = 0.0001;

    m_ReferenceImageResampler = 0;
    m_MovingImageResampler = 0;

    m_VerboseProgression = true;

    this->SetNumberOfThreads(this->GetMultiThreader()->GetNumberOfThreads());

    m_InitialTransform = 0;
    this->SetNumberOfRequiredOutputs(1);
    TransformOutputPointer transformDecorator = static_cast <TransformOutputType *> (this->MakeOutput(0).GetPointer());
    this->itk::ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());
}

/**
 *  Get Output
 */
template <typename TInputImageType>
typename BaseBMRegistrationMethod <TInputImageType>::TransformOutputType *
BaseBMRegistrationMethod <TInputImageType>
::GetOutput()
{
    return static_cast <TransformOutputType *> (this->ProcessObject::GetOutput(0));
}

template <typename TInputImageType>
itk::DataObject::Pointer
BaseBMRegistrationMethod <TInputImageType>
::MakeOutput(DataObjectPointerArraySizeType output)
{
    switch (output)
    {
        case 0:
            return static_cast <itk::DataObject*> (TransformOutputType::New().GetPointer());
            break;
        default:
            itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
            return 0;
    }
}

/**
 * Starts registration when Update method is called
 */
template <typename TInputImageType>
void
BaseBMRegistrationMethod <TInputImageType>
::GenerateData()
{
    m_Abort = false;
    this->StartOptimization();
}

/**
 * Starts the Registration Process
 */
template <typename TInputImageType>
void
BaseBMRegistrationMethod <TInputImageType>
::StartOptimization()
{
    m_Agregator->SetInputTransformType(m_BlockMatcher->GetAgregatorInputTransformType());

    TransformPointer computedTransform = NULL;
    this->SetupTransform(computedTransform);

    //progress management
    itk::ProgressReporter progress(this, 0, m_MaximumIterations);

    // Real work goes here
    InputImagePointer fixedResampled, movingResampled;
    for (unsigned int iterations = 0; iterations < m_MaximumIterations && !m_Abort; ++iterations)
    {
        // Resample fixed and moving image here
        this->ResampleImages(computedTransform, fixedResampled, movingResampled);

        // Perform one iteration of registration between the images
        // Calls pure virtual method that can use the block matching class available here
        this->GetAgregator()->SetCurrentLinearTransform(computedTransform);
        TransformPointer addOn;
        this->PerformOneIteration(fixedResampled, movingResampled, addOn);

        bool continueLoop = this->ComposeAddOnWithTransform(computedTransform,addOn);

        if (m_VerboseProgression)
            std::cout << "Iteration " << iterations << " done..." << std::endl;

        if (iterations != m_MaximumIterations - 1)
            progress.CompletedPixel();

        if (!continueLoop)
            break;
    }

    TransformOutputPointer transformDecorator = TransformOutputType::New();
    transformDecorator->Set(computedTransform.GetPointer());

    this->itk::ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());
}

template <typename TInputImageType>
void
BaseBMRegistrationMethod <TInputImageType>
::SetupTransform(TransformPointer &optimizedTransform)
{
    if (m_Agregator->GetOutputTransformType() != AgregatorType::SVF)
    {
        if (m_InitialTransform)
        {
            optimizedTransform = AffineTransformType::New();
            optimizedTransform->SetParameters(m_InitialTransform->GetParameters());
        }
        else
        {
            typename AffineTransformType::Pointer tmpTrsf = AffineTransformType::New();
            tmpTrsf->SetIdentity();
            optimizedTransform = tmpTrsf;
        }
    }
    else
    {
        if (m_InitialTransform)
            optimizedTransform = m_InitialTransform;
        else
        {
            SVFTransformPointer tmpTrsf = SVFTransformType::New();
            tmpTrsf->SetIdentity();
            optimizedTransform = tmpTrsf;
        }
    }
}

template <typename TInputImageType>
void
BaseBMRegistrationMethod <TInputImageType>
::ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &movingImage)
{
    // Moving image resampling
    if (m_MovingImage->GetNumberOfComponentsPerPixel() > 1)
    {
        // Model resampling
        typedef anima::OrientedModelBaseResampleImageFilter <InputImageType,AgregatorScalarType> InternalFilterType;
        InternalFilterType *resampleFilter = dynamic_cast <InternalFilterType *> (m_MovingImageResampler.GetPointer());

        if (m_Agregator->GetOutputTransformType() == AgregatorType::SVF)
        {
            // Compute temporary field and set it to resampler
            DisplacementFieldTransformPointer dispTrsf = DisplacementFieldTransformType::New();
            SVFTransformType *svfCast = dynamic_cast<SVFTransformType *> (currentTransform);

            anima::GetSVFExponential(svfCast,dispTrsf.GetPointer(),m_ExponentiationOrder,this->GetNumberOfThreads(),false);

            resampleFilter->SetTransform(dispTrsf);
        }
        else
            resampleFilter->SetTransform(currentTransform);
    }
    else
    {
        // Anatomical resampling
        typedef itk::Image <ImageScalarType, TInputImageType::ImageDimension> InternalScalarImageType;
        typedef anima::ResampleImageFilter <InternalScalarImageType,InternalScalarImageType,AgregatorScalarType> InternalFilterType;
        InternalFilterType *resampleFilter = dynamic_cast <InternalFilterType *> (m_MovingImageResampler.GetPointer());

        if (m_Agregator->GetOutputTransformType() == AgregatorType::SVF)
        {
            // Compute temporary field and set it to resampler
            DisplacementFieldTransformPointer dispTrsf = DisplacementFieldTransformType::New();
            SVFTransformType *svfCast = dynamic_cast<SVFTransformType *> (currentTransform);

            anima::GetSVFExponential(svfCast,dispTrsf.GetPointer(),m_ExponentiationOrder,this->GetNumberOfThreads(),false);

            resampleFilter->SetTransform(dispTrsf);
        }
        else
            resampleFilter->SetTransform(currentTransform);
    }

    m_MovingImageResampler->SetInput(m_MovingImage);
    m_MovingImageResampler->SetNumberOfThreads(this->GetNumberOfThreads());
    m_MovingImageResampler->Update();

    movingImage = m_MovingImageResampler->GetOutput();
    movingImage->DisconnectPipeline();

    // Fixed image resampling
    if (m_FixedImage->GetNumberOfComponentsPerPixel() > 1)
    {
        // Model resampling
        typedef anima::OrientedModelBaseResampleImageFilter <InputImageType,AgregatorScalarType> InternalFilterType;
        InternalFilterType *resampleFilter = dynamic_cast <InternalFilterType *> (m_ReferenceImageResampler.GetPointer());

        if (m_Agregator->GetOutputTransformType() == AgregatorType::SVF)
        {
            // Compute temporary field and set it to resampler
            DisplacementFieldTransformPointer dispTrsf = DisplacementFieldTransformType::New();
            SVFTransformType *svfCast = dynamic_cast<SVFTransformType *> (currentTransform);

            anima::GetSVFExponential(svfCast,dispTrsf.GetPointer(),m_ExponentiationOrder,this->GetNumberOfThreads(),true);

            resampleFilter->SetTransform(dispTrsf);
        }
        else
        {
            AffineTransformType *affCast = dynamic_cast<AffineTransformType *> (currentTransform);
            TransformPointer reverseTrsf = affCast->GetInverseTransform();
            resampleFilter->SetTransform(reverseTrsf);
        }
    }
    else
    {
        // Anatomical resampling
        typedef itk::Image <ImageScalarType, TInputImageType::ImageDimension> InternalScalarImageType;
        typedef anima::ResampleImageFilter <InternalScalarImageType,InternalScalarImageType,AgregatorScalarType> InternalFilterType;
        InternalFilterType *resampleFilter = dynamic_cast <InternalFilterType *> (m_ReferenceImageResampler.GetPointer());

        if (m_Agregator->GetOutputTransformType() == AgregatorType::SVF)
        {
            // Compute temporary field and set it to resampler
            DisplacementFieldTransformPointer dispTrsf = DisplacementFieldTransformType::New();
            SVFTransformType *svfCast = dynamic_cast<SVFTransformType *> (currentTransform);

            anima::GetSVFExponential(svfCast,dispTrsf.GetPointer(),m_ExponentiationOrder,this->GetNumberOfThreads(),true);

            resampleFilter->SetTransform(dispTrsf);
        }
        else
        {
            AffineTransformType *affCast = dynamic_cast<AffineTransformType *> (currentTransform);
            TransformPointer reverseTrsf = affCast->GetInverseTransform();
            resampleFilter->SetTransform(reverseTrsf);
        }
    }

    m_ReferenceImageResampler->SetInput(m_FixedImage);
    m_ReferenceImageResampler->SetNumberOfThreads(this->GetNumberOfThreads());
    m_ReferenceImageResampler->Update();

    refImage = m_ReferenceImageResampler->GetOutput();
    refImage->DisconnectPipeline();
}

template <typename TInputImageType>
bool
BaseBMRegistrationMethod <TInputImageType>
::ComposeAddOnWithTransform(TransformPointer &computedTransform, TransformType *addOn)
{
    if (m_Agregator->GetOutputTransformType() != AgregatorType::SVF)
    {
        typename TransformType::ParametersType oldPars = computedTransform->GetParameters();

        if (m_Agregator->GetOutputTransformType() != AgregatorType::ANISOTROPIC_SIM)
        {
            AffineTransformType *tmpTrsf = dynamic_cast<AffineTransformType *>(computedTransform.GetPointer());
            AffineTransformType *tmpAddOn = dynamic_cast<AffineTransformType *>(addOn);
            tmpTrsf->Compose(tmpAddOn, true);
        }
        else
        {
            // no composition since we compute global transform only
            computedTransform->SetParameters(addOn->GetParameters());
        }

        typename TransformType::ParametersType newPars = computedTransform->GetParameters();

        // Compute the distance between consecutive solutions, until a certain threshold
        double err = 0;
        for (unsigned int i = 0; i < newPars.Size(); ++i)
            err += std::pow(newPars[i] - oldPars[i], 2.);

        if (err <= m_MinimalTransformError)
            return false;
    }
    else
    {
        // Add update to current velocity field (cf. Vercauteren et al, 2008)
        SVFTransformType *tmpTrsf = dynamic_cast<SVFTransformType *>(computedTransform.GetPointer());
        SVFTransformType *tmpAddOn = dynamic_cast<SVFTransformType *>(addOn);

        anima::composeSVF(tmpTrsf,tmpAddOn,this->GetNumberOfThreads(),m_BCHCompositionOrder);

        typedef typename SVFTransformType::VectorFieldType VectorFieldType;
        typedef itk::ImageRegionConstIterator <VectorFieldType> IteratorType;
        IteratorType diffItr(tmpAddOn->GetParametersAsVectorField(),
                             tmpAddOn->GetParametersAsVectorField()->GetLargestPossibleRegion());

        bool smallEnoughTransform = true;
        while (!diffItr.IsAtEnd())
        {
            double err = 0;
            for (unsigned int i = 0;i < VectorFieldType::ImageDimension;++i)
                err += diffItr.Get()[i] * diffItr.Get()[i];

            if (err > m_MinimalTransformError)
            {
                smallEnoughTransform = false;
                break;
            }

            ++diffItr;
        }

        if (smallEnoughTransform)
            return false;

        if (m_SVFElasticRegSigma > 0)
        {
            typedef typename SVFTransformType::VectorFieldType VelocityFieldType;

            typedef anima::SmoothingRecursiveYvvGaussianImageFilter <VelocityFieldType, VelocityFieldType> SmoothingFilterType;
            typename SmoothingFilterType::Pointer smootherPtr = SmoothingFilterType::New();

            smootherPtr->SetInput(tmpTrsf->GetParametersAsVectorField());
            smootherPtr->SetSigma(m_SVFElasticRegSigma);
            smootherPtr->SetNumberOfThreads(this->GetNumberOfThreads());

            smootherPtr->Update();

            typename VelocityFieldType::Pointer tmpSmoothed = smootherPtr->GetOutput();
            tmpSmoothed->DisconnectPipeline();

            tmpTrsf->SetParametersAsVectorField(tmpSmoothed);
        }
    }

    return true;
}

/**
 * PrintSelf
 */
template <typename TInputImageType>
void
BaseBMRegistrationMethod <TInputImageType>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf( os, indent );

    os << indent << "Fixed Image: " << m_FixedImage.GetPointer() << std::endl;
    os << indent << "Moving Image: " << m_MovingImage.GetPointer() << std::endl;

    os << indent << "Maximum Iterations: " << m_MaximumIterations << std::endl;
}

} // end namespace anima
