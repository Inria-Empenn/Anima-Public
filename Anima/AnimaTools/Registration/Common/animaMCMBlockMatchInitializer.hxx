#pragma once
#include "animaMCMBlockMatchInitializer.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkExpNegativeImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>

namespace anima
{

template <class PixelType, unsigned int NDimensions>
void
MCMBlockMatchingInitializer<PixelType,NDimensions>
::AddReferenceImage(itk::ImageBase <NDimensions> *refImage)
{
    this->Superclass::AddReferenceImage(refImage);

    MCMImageType *castImage = dynamic_cast <MCMImageType *> (refImage);
    if (castImage != 0)
        m_ReferenceModels.push_back(castImage->GetDescriptionModel());
}

template <class PixelType, unsigned int NDimensions>
void
MCMBlockMatchingInitializer<PixelType,NDimensions>
::InitializeThreading(unsigned int maskIndex, BlockGeneratorThreadStruct *&workStr)
{
    MCMBlockGeneratorThreadStruct *tmpStr = new MCMBlockGeneratorThreadStruct;
    workStr = tmpStr;
    this->Superclass::InitializeThreading(maskIndex,workStr);

    tmpStr->ref_models.resize(this->GetNumberOfThreads());
    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        std::vector <MCModelPointer> ref_models(m_ReferenceModels.size());
        for (unsigned int j = 0;j < m_ReferenceModels.size();++j)
            ref_models[j] = m_ReferenceModels[j]->Clone();

        tmpStr->ref_models[i] = ref_models;
    }
}

template <class PixelType, unsigned int NDimensions>
bool MCMBlockMatchingInitializer<PixelType,NDimensions>
::CheckOrientedModelVariance(unsigned int imageIndex, ImageRegionType &region, double &blockVariance,
                             BlockGeneratorThreadStruct *workStr, unsigned int threadId)
{
    MCMImageType *refImage = dynamic_cast <MCMImageType *> (this->GetReferenceVectorImage(imageIndex));
    MCMBlockGeneratorThreadStruct *tmpStr = dynamic_cast <MCMBlockGeneratorThreadStruct *> (workStr);

    itk::ImageRegionConstIterator <MCMImageType> refItr(refImage,region);
    typedef typename MCMImageType::PixelType VectorType;

    MCModelPointer refModel = tmpStr->ref_models[threadId][imageIndex];
    unsigned int nbPts = 0;
    unsigned int numIsoCompartments = refModel->GetNumberOfIsotropicCompartments();

    double meanVal = 0;
    blockVariance = 0;
    while (!refItr.IsAtEnd())
    {
        refModel->SetModelVector(refItr.Get());

        // Here we sum all isotropic compartment weights
        double tmpVal = 0;
        for (unsigned int i = 0;i < numIsoCompartments;++i)
            tmpVal += refModel->GetCompartmentWeight(i);

        meanVal += tmpVal;
        blockVariance += tmpVal * tmpVal;

        ++nbPts;
        ++refItr;
    }

    if (nbPts <= 1)
        return false;

    blockVariance = (blockVariance - meanVal * meanVal / nbPts) / (nbPts - 1.0);

    if (blockVariance > this->GetOrientedModelVarianceThreshold())
        return true;
    else
        return false;
}

}// end of namespace anima
