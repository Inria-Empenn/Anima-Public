#pragma once
#include "animaVelocityUtils.h"

#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>

#include <itkComposeDisplacementFieldsImageFilter.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <animaSVFLieBracketImageFilter.h>
#include <animaSVFExponentialImageFilter.h>

namespace anima
{

template <class ScalarType, unsigned int NDimensions>
void composeSVF(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *addonTrsf,
                unsigned int numThreads, unsigned int bchOrder)
{
    if ((baseTrsf->GetParametersAsVectorField() == NULL)&&(addonTrsf->GetParametersAsVectorField() == NULL))
        return;

    if ((bchOrder > 4)||(bchOrder < 1))
        throw itk::ExceptionObject(__FILE__,__LINE__,"Invalid BCH order, not implemented yet",ITK_LOCATION);

    typedef typename itk::StationaryVelocityFieldTransform <ScalarType,NDimensions>::VectorFieldType VelocityFieldType;

    if (baseTrsf->GetParametersAsVectorField() == NULL)
    {
        baseTrsf->SetParametersAsVectorField(addonTrsf->GetParametersAsVectorField());
        return;
    }

    typedef itk::AddImageFilter <VelocityFieldType, VelocityFieldType> AddFilterType;    
    typedef itk::MultiplyImageFilter <VelocityFieldType, itk::Image <double, NDimensions>,
            VelocityFieldType> MultiplyConstFilterType;

    typename AddFilterType::Pointer bchAdder = AddFilterType::New();
    bchAdder->SetInput(0,baseTrsf->GetParametersAsVectorField());
    bchAdder->SetInput(1,addonTrsf->GetParametersAsVectorField());

    if (numThreads > 0)
        bchAdder->SetNumberOfWorkUnits(numThreads);

    bchAdder->Update();

    typename VelocityFieldType::Pointer resField = bchAdder->GetOutput();
    resField->DisconnectPipeline();

    typedef anima::SVFLieBracketImageFilter <ScalarType, NDimensions> LieBracketFilterType;
    typename LieBracketFilterType::JacobianImagePointer baseTrsfJac, addonTrsfJac;
    typename LieBracketFilterType::OutputImagePointer previousLieBracket;

    if (bchOrder >= 2)
    {
        // Compute Lie bracket and add half of it to the output
        typename LieBracketFilterType::Pointer lieBracketFilter = LieBracketFilterType::New();

        lieBracketFilter->SetInput(0,baseTrsf->GetParametersAsVectorField());
        lieBracketFilter->SetInput(1,addonTrsf->GetParametersAsVectorField());

        if (numThreads > 0)
            lieBracketFilter->SetNumberOfWorkUnits(numThreads);

        lieBracketFilter->Update();
        baseTrsfJac = lieBracketFilter->GetFirstFieldJacobian();
        baseTrsfJac->DisconnectPipeline();
        addonTrsfJac = lieBracketFilter->GetSecondFieldJacobian();
        addonTrsfJac->DisconnectPipeline();
        previousLieBracket = lieBracketFilter->GetOutput();
        previousLieBracket->DisconnectPipeline();

        typename MultiplyConstFilterType::Pointer bracketMultiplier = MultiplyConstFilterType::New();
        bracketMultiplier->SetInput(0,previousLieBracket);
        bracketMultiplier->SetConstant(0.5);

        if (numThreads > 0)
            bracketMultiplier->SetNumberOfWorkUnits(numThreads);

        bracketMultiplier->Update();

        typename AddFilterType::Pointer bchSecondAdder = AddFilterType::New();
        bchSecondAdder->SetInput(0,resField);
        bchSecondAdder->SetInput(1,bracketMultiplier->GetOutput());

        if (numThreads > 0)
            bchSecondAdder->SetNumberOfWorkUnits(numThreads);

        bchSecondAdder->Update();

        resField = bchSecondAdder->GetOutput();
        resField->DisconnectPipeline();
    }

    if (bchOrder >= 3)
    {
        // Compute Lie bracket one way and add 1/12 of it to the output
        typename LieBracketFilterType::Pointer lieBracketFilter = LieBracketFilterType::New();

        lieBracketFilter->SetInput(0,baseTrsf->GetParametersAsVectorField());
        lieBracketFilter->SetInput(1,previousLieBracket);
        lieBracketFilter->SetFirstFieldJacobian(baseTrsfJac);

        if (numThreads > 0)
            lieBracketFilter->SetNumberOfWorkUnits(numThreads);

        lieBracketFilter->Update();

        typename MultiplyConstFilterType::Pointer bracketMultiplier = MultiplyConstFilterType::New();
        bracketMultiplier->SetInput(0,lieBracketFilter->GetOutput());
        bracketMultiplier->SetConstant(1.0 / 12);

        if (numThreads > 0)
            bracketMultiplier->SetNumberOfWorkUnits(numThreads);

        bracketMultiplier->Update();

        typename AddFilterType::Pointer bchSecondAdder = AddFilterType::New();
        bchSecondAdder->SetInput(0,resField);
        bchSecondAdder->SetInput(1,bracketMultiplier->GetOutput());

        if (numThreads > 0)
            bchSecondAdder->SetNumberOfWorkUnits(numThreads);

        bchSecondAdder->Update();

        resField = bchSecondAdder->GetOutput();
        resField->DisconnectPipeline();

        // Compute Lie bracket the other way round and add 1/12 of it to the output
        typename LieBracketFilterType::Pointer reverseLieBracketFilter = LieBracketFilterType::New();

        reverseLieBracketFilter->SetInput(0,previousLieBracket);
        reverseLieBracketFilter->SetInput(1,addonTrsf->GetParametersAsVectorField());
        reverseLieBracketFilter->SetFirstFieldJacobian(lieBracketFilter->GetSecondFieldJacobian());
        reverseLieBracketFilter->SetSecondFieldJacobian(addonTrsfJac);

        if (numThreads > 0)
            reverseLieBracketFilter->SetNumberOfWorkUnits(numThreads);

        reverseLieBracketFilter->Update();

        bracketMultiplier = MultiplyConstFilterType::New();
        bracketMultiplier->SetInput(0,reverseLieBracketFilter->GetOutput());
        bracketMultiplier->SetConstant(1.0 / 12);

        if (numThreads > 0)
            bracketMultiplier->SetNumberOfWorkUnits(numThreads);

        bracketMultiplier->Update();

        bchSecondAdder = AddFilterType::New();
        bchSecondAdder->SetInput(0,resField);
        bchSecondAdder->SetInput(1,bracketMultiplier->GetOutput());

        if (numThreads > 0)
            bchSecondAdder->SetNumberOfWorkUnits(numThreads);

        bchSecondAdder->Update();

        resField = bchSecondAdder->GetOutput();
        resField->DisconnectPipeline();

        previousLieBracket = lieBracketFilter->GetOutput();
        previousLieBracket->DisconnectPipeline();
    }

    if (bchOrder == 4)
    {
        // Compute Lie bracket one way and add 1/24 of it to the output
        typename LieBracketFilterType::Pointer lieBracketFilter = LieBracketFilterType::New();

        lieBracketFilter->SetInput(0,previousLieBracket);
        lieBracketFilter->SetInput(1,addonTrsf->GetParametersAsVectorField());
        lieBracketFilter->SetSecondFieldJacobian(addonTrsfJac);

        if (numThreads > 0)
            lieBracketFilter->SetNumberOfWorkUnits(numThreads);

        lieBracketFilter->Update();

        typename MultiplyConstFilterType::Pointer bracketMultiplier = MultiplyConstFilterType::New();
        bracketMultiplier->SetInput(0,lieBracketFilter->GetOutput());
        bracketMultiplier->SetConstant(1.0 / 24);

        if (numThreads > 0)
            bracketMultiplier->SetNumberOfWorkUnits(numThreads);

        bracketMultiplier->Update();

        typename AddFilterType::Pointer bchSecondAdder = AddFilterType::New();
        bchSecondAdder->SetInput(0,resField);
        bchSecondAdder->SetInput(1,bracketMultiplier->GetOutput());

        if (numThreads > 0)
            bchSecondAdder->SetNumberOfWorkUnits(numThreads);

        bchSecondAdder->Update();

        resField = bchSecondAdder->GetOutput();
        resField->DisconnectPipeline();
    }

    baseTrsf->SetParametersAsVectorField(resField.GetPointer());
}

template <class ScalarType, unsigned int NDimensions>
void GetSVFExponential(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                       rpi::DisplacementFieldTransform <ScalarType,NDimensions> *resultTransform,
                       unsigned int exponentiationOrder, unsigned int numThreads, bool invert)
{
    if (baseTrsf->GetParametersAsVectorField() == NULL)
        return;

    typedef itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> SVFType;
    typedef typename SVFType::VectorFieldType FieldType;
    typedef typename FieldType::Pointer FieldPointer;

    typedef anima::SVFExponentialImageFilter <ScalarType, NDimensions> ExponentialFilterType;

    FieldPointer tmpPtr = const_cast <FieldType *> (baseTrsf->GetParametersAsVectorField());
    if (invert)
    {
        typedef itk::MultiplyImageFilter <FieldType,itk::Image <double, NDimensions>, FieldType> MultiplyFilterType;
        typename MultiplyFilterType::Pointer multiplier = MultiplyFilterType::New();
        multiplier->SetInput(tmpPtr);
        multiplier->SetConstant(-1.0);

        multiplier->SetNumberOfWorkUnits(numThreads);
        multiplier->Update();

        tmpPtr = multiplier->GetOutput();
        tmpPtr->DisconnectPipeline();
    }

    typename ExponentialFilterType::Pointer expFilter = ExponentialFilterType::New();
    expFilter->SetInput(tmpPtr);
    expFilter->SetExponentiationOrder(exponentiationOrder);
    expFilter->SetNumberOfWorkUnits(numThreads);
    expFilter->SetMaximalDisplacementAmplitude(0.25);

    expFilter->Update();

    FieldPointer resField = expFilter->GetOutput();
    resField->DisconnectPipeline();

    resultTransform->SetParametersAsVectorField(resField.GetPointer());
}

template <class ScalarType, unsigned int NDimensions>
void composeDistortionCorrections(typename rpi::DisplacementFieldTransform <ScalarType,NDimensions>::Pointer &baseTrsf,
                                  typename rpi::DisplacementFieldTransform <ScalarType,NDimensions>::Pointer &positiveAddOn,
                                  typename rpi::DisplacementFieldTransform <ScalarType,NDimensions>::Pointer &negativeAddOn,
                                  unsigned int numThreads)
{
    typedef rpi::DisplacementFieldTransform <ScalarType,NDimensions> DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::VectorFieldType VectorFieldType;
    typedef itk::ComposeDisplacementFieldsImageFilter <VectorFieldType,VectorFieldType> ComposeFilterType;
    typedef itk::MultiplyImageFilter <VectorFieldType,itk::Image <double, NDimensions>, VectorFieldType> MultiplyFilterType;
    typedef typename itk::ImageRegionIterator <VectorFieldType> VectorFieldIterator;
    typedef typename VectorFieldType::PixelType VectorType;

    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction <VectorFieldType,
            typename VectorType::ValueType> VectorInterpolateFunctionType;

    typename ComposeFilterType::Pointer composePositiveFilter = ComposeFilterType::New();
    composePositiveFilter->SetWarpingField(positiveAddOn->GetParametersAsVectorField());
    composePositiveFilter->SetDisplacementField(baseTrsf->GetParametersAsVectorField());
    composePositiveFilter->SetNumberOfWorkUnits(numThreads);

    typename VectorInterpolateFunctionType::Pointer interpolator = VectorInterpolateFunctionType::New();

    composePositiveFilter->SetInterpolator(interpolator);
    composePositiveFilter->Update();
    positiveAddOn->SetParametersAsVectorField(composePositiveFilter->GetOutput());

    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    multiplyFilter->SetInput(baseTrsf->GetParametersAsVectorField());
    multiplyFilter->SetConstant(-1.0);

    if (numThreads > 0)
        multiplyFilter->SetNumberOfWorkUnits(numThreads);

    multiplyFilter->InPlaceOn();

    multiplyFilter->Update();

    typename ComposeFilterType::Pointer composeNegativeFilter = ComposeFilterType::New();
    composeNegativeFilter->SetWarpingField(negativeAddOn->GetParametersAsVectorField());
    composeNegativeFilter->SetDisplacementField(multiplyFilter->GetOutput());

    if (numThreads > 0)
        composeNegativeFilter->SetNumberOfWorkUnits(numThreads);

    interpolator = VectorInterpolateFunctionType::New();

    composeNegativeFilter->SetInterpolator(interpolator);
    composeNegativeFilter->Update();
    negativeAddOn->SetParametersAsVectorField(composeNegativeFilter->GetOutput());

    VectorFieldIterator positiveItr(const_cast <VectorFieldType *> (positiveAddOn->GetParametersAsVectorField()),
                                    positiveAddOn->GetParametersAsVectorField()->GetLargestPossibleRegion());

    VectorFieldIterator negativeItr(const_cast <VectorFieldType *> (negativeAddOn->GetParametersAsVectorField()),
                                    negativeAddOn->GetParametersAsVectorField()->GetLargestPossibleRegion());

    // And compose them to get a transformation conform to disto correction requirements
    VectorType tmpVec;
    while (!positiveItr.IsAtEnd())
    {
        tmpVec = 0.5 * (positiveItr.Get() - negativeItr.Get());
        positiveItr.Set(tmpVec);
        negativeItr.Set(- tmpVec);

        ++positiveItr;
        ++negativeItr;
    }

    baseTrsf = positiveAddOn;
}

} // end of namespace anima
