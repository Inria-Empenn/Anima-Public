#pragma once
#include "animaVelocityUtils.h"

#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>

#include <itkComposeDisplacementFieldsImageFilter.h>

namespace anima
{

template <class ScalarType, unsigned int NDimensions>
void composeSVF(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *addonTrsf,
                unsigned int numThreads)
{
    if ((baseTrsf->GetParametersAsVectorField() == NULL)&&(addonTrsf->GetParametersAsVectorField() == NULL))
        return;

    typedef typename itk::StationaryVelocityFieldTransform <ScalarType,NDimensions>::VectorFieldType VelocityFieldType;

    if (baseTrsf->GetParametersAsVectorField() == NULL)
    {
        baseTrsf->SetParametersAsVectorField(addonTrsf->GetParametersAsVectorField());
        return;
    }

    typedef itk::AddImageFilter <VelocityFieldType, VelocityFieldType> AddFilterType;
    typename AddFilterType::Pointer bchAdder = AddFilterType::New();
    bchAdder->SetInput(0,baseTrsf->GetParametersAsVectorField());
    bchAdder->SetInput(1,addonTrsf->GetParametersAsVectorField());

    if (numThreads > 0)
        bchAdder->SetNumberOfThreads(numThreads);

    bchAdder->Update();

    typename VelocityFieldType::Pointer resField = bchAdder->GetOutput();
    resField->DisconnectPipeline();

    baseTrsf->SetParametersAsVectorField(resField.GetPointer());
}

template <class ScalarType, unsigned int NDimensions>
void GetSVFExponential(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                       rpi::DisplacementFieldTransform <ScalarType,NDimensions> *resultTransform, bool invert)
{
    if (baseTrsf->GetParametersAsVectorField() == NULL)
        return;

    typedef itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> SVFType;
    typedef typename SVFType::VectorFieldType FieldType;
    typedef typename FieldType::Pointer FieldPointer;

    typename SVFType::Pointer tmpPtr = baseTrsf;
    if (invert)
        tmpPtr = dynamic_cast <SVFType *> (baseTrsf->GetInverseTransform().GetPointer());

    FieldPointer resField = tmpPtr->GetDisplacementFieldAsVectorField();
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
    typedef itk::MultiplyImageFilter <VectorFieldType,itk::Image <float, NDimensions>, VectorFieldType> MultiplyFilterType;
    typedef typename itk::ImageRegionIterator <VectorFieldType> VectorFieldIterator;
    typedef typename VectorFieldType::PixelType VectorType;

    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction <VectorFieldType,
            typename VectorType::ValueType> VectorInterpolateFunctionType;

    typename ComposeFilterType::Pointer composePositiveFilter = ComposeFilterType::New();
    composePositiveFilter->SetWarpingField(positiveAddOn->GetParametersAsVectorField());
    composePositiveFilter->SetDisplacementField(baseTrsf->GetParametersAsVectorField());
    composePositiveFilter->SetNumberOfThreads(numThreads);

    typename VectorInterpolateFunctionType::Pointer interpolator = VectorInterpolateFunctionType::New();

    composePositiveFilter->SetInterpolator(interpolator);
    composePositiveFilter->Update();
    positiveAddOn->SetParametersAsVectorField(composePositiveFilter->GetOutput());

    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    multiplyFilter->SetInput(baseTrsf->GetParametersAsVectorField());
    multiplyFilter->SetConstant(-1.0);

    if (numThreads > 0)
        multiplyFilter->SetNumberOfThreads(numThreads);

    multiplyFilter->InPlaceOn();

    multiplyFilter->Update();

    typename ComposeFilterType::Pointer composeNegativeFilter = ComposeFilterType::New();
    composeNegativeFilter->SetWarpingField(negativeAddOn->GetParametersAsVectorField());
    composeNegativeFilter->SetDisplacementField(multiplyFilter->GetOutput());

    if (numThreads > 0)
        composeNegativeFilter->SetNumberOfThreads(numThreads);

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
