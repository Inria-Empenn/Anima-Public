#pragma once
#include "animaBalooSVFTransformAgregator.h"

#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>

#include <animaLogRigid3DTransform.h>
#include <animaDirectionScaleSkewTransform.h>

#include <animaMatrixLogExp.h>
#include <itkMultiThreader.h>
#include <itkTimeProbe.h>
#include <itkImageFileWriter.h>

#include <boost/math/special_functions/factorials.hpp>

namespace anima
{

template <unsigned int NDimensions>
BalooSVFTransformAgregator <NDimensions>::
BalooSVFTransformAgregator() : Superclass()
{
    m_ExtrapolationSigma = 4.0;
    m_OutlierRejectionSigma = 3.0;

    m_NumberOfThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
}

template <unsigned int NDimensions>
bool
BalooSVFTransformAgregator <NDimensions>::
Update()
{
    this->SetUpToDate(false);

    if (this->GetInputWeights().size() != this->GetInputTransforms().size())
        return false;

    if (this->GetOutputTransformType() != Superclass::SVF)
        return false;

    itk::TimeProbe tmpTime;
    tmpTime.Start();
    switch (this->GetInputTransformType())
    {
        case Superclass::TRANSLATION:
            this->estimateSVFFromTranslations();
            this->SetUpToDate(true);
            break;

        case Superclass::RIGID:
            this->estimateSVFFromRigidTransforms();
            this->SetUpToDate(true);
            break;

        case Superclass::AFFINE:
        default:
            this->estimateSVFFromAffineTransforms();
            this->SetUpToDate(true);
            break;
    }

    tmpTime.Stop();
    std::cout << "Agregation performed in " << tmpTime.GetTotal() << std::endl;

    return true;
}

template <unsigned int NDimensions>
void
BalooSVFTransformAgregator <NDimensions>::
estimateSVFFromTranslations()
{
    unsigned int nbPts = this->GetInputRegions().size();
    VelocityFieldPointer velocityField = VelocityFieldType::New();
    velocityField->Initialize();

    velocityField->SetRegions (m_LargestRegion);
    velocityField->SetSpacing (m_Spacing);
    velocityField->SetOrigin (m_Origin);
    velocityField->SetDirection (m_Direction);
    velocityField->Allocate();

    VelocityFieldPixelType zeroDisp;
    zeroDisp.Fill(0);
    itk::ImageRegionIterator < VelocityFieldType > svfIterator(velocityField,m_LargestRegion);
    while (!svfIterator.IsAtEnd())
    {
        svfIterator.Set(zeroDisp);
        ++svfIterator;
    }

    WeightImagePointer weights = WeightImageType::New();
    weights->Initialize();

    weights->SetRegions (m_LargestRegion);
    weights->SetSpacing (m_Spacing);
    weights->SetOrigin (m_Origin);
    weights->SetDirection (m_Direction);
    weights->Allocate();
    weights->FillBuffer(anima::zeroWeight);

    ParametersType tmpParams;
    std::vector <VelocityFieldPixelType> curDisps(nbPts);
    std::vector <VelocityFieldIndexType> posIndexes(nbPts);

    for (unsigned int i = 0;i < nbPts;++i)
    {
        tmpParams = this->GetInputTransform(i)->GetParameters();
        double tmpWeight = this->GetInputWeight(i);
        weights->TransformPhysicalPointToIndex(this->GetInputOrigin(i),posIndexes[i]);

        for (unsigned int j = 0;j < NDimensions;++j)
            curDisps[i][j] = tmpParams[j];

        velocityField->SetPixel(posIndexes[i],curDisps[i] * tmpWeight);
        weights->SetPixel(posIndexes[i],tmpWeight);
    }

    for (unsigned int i = 0;i < m_DamIndexes.size();++i)
        weights->SetPixel(m_DamIndexes[i],1.0);

    anima::filterInputs<ScalarType,NDimensions,NDimensions>(weights,velocityField,curDisps,posIndexes,this);

    // Create the final transform
    typename BaseOutputTransformType::Pointer resultTransform = BaseOutputTransformType::New();
    resultTransform->SetIdentity();
    resultTransform->SetParametersAsVectorField(velocityField.GetPointer());

    this->SetOutput(resultTransform);
}

template <unsigned int NDimensions>
void
BalooSVFTransformAgregator <NDimensions>::
estimateSVFFromRigidTransforms()
{
    const unsigned int NDegreesFreedom = NDimensions * (NDimensions + 1) / 2;
    typedef itk::Image < itk::Vector <ScalarType, NDegreesFreedom>, NDimensions > RigidFieldType;
    typedef typename RigidFieldType::Pointer RigidFieldPointer;
    typedef itk::Vector <ScalarType, NDegreesFreedom> RigidVectorType;

    unsigned int nbPts = this->GetInputRegions().size();
    RigidFieldPointer rigidField = RigidFieldType::New();
    rigidField->Initialize();

    rigidField->SetRegions (m_LargestRegion);
    rigidField->SetSpacing (m_Spacing);
    rigidField->SetOrigin (m_Origin);
    rigidField->SetDirection (m_Direction);
    rigidField->Allocate();

    RigidVectorType zeroDisp;
    zeroDisp.Fill(0);
    itk::ImageRegionIterator < RigidFieldType > rigidIterator(rigidField,m_LargestRegion);
    while (!rigidIterator.IsAtEnd())
    {
        rigidIterator.Set(zeroDisp);
        ++rigidIterator;
    }

    WeightImagePointer weights = WeightImageType::New();
    weights->Initialize();

    weights->SetRegions (m_LargestRegion);
    weights->SetSpacing (m_Spacing);
    weights->SetOrigin (m_Origin);
    weights->SetDirection (m_Direction);
    weights->Allocate();
    weights->FillBuffer(anima::zeroWeight);

    RigidVectorType curLog;

    vnl_matrix <InternalScalarType> tmpMatrix(NDimensions+1,NDimensions+1,0);
    tmpMatrix(NDimensions,NDimensions) = 1;
    std::vector < RigidVectorType > logVectors(nbPts);
    std::vector <VelocityFieldIndexType> posIndexes(nbPts);
    VelocityFieldPointType curPoint;

    typedef anima::LogRigid3DTransform <InternalScalarType> LogRigidTransformType;

    for (unsigned int i = 0;i < nbPts;++i)
    {
        double tmpWeight = this->GetInputWeight(i);
        weights->TransformPhysicalPointToIndex(this->GetInputOrigin(i),posIndexes[i]);

        LogRigidTransformType *tmpTrsf = (LogRigidTransformType *)this->GetInputTransform(i);
        logVectors[i] = tmpTrsf->GetLogVector();
        rigidField->SetPixel(posIndexes[i],logVectors[i] * tmpWeight);
        weights->SetPixel(posIndexes[i],tmpWeight);
    }

    for (unsigned int i = 0;i < m_DamIndexes.size();++i)
        weights->SetPixel(m_DamIndexes[i],1.0);

    anima::filterInputs<ScalarType,NDimensions,NDegreesFreedom>(weights,rigidField,logVectors,posIndexes,this);

    VelocityFieldPointer velocityField = VelocityFieldType::New();
    velocityField->Initialize();

    velocityField->SetRegions (m_LargestRegion);
    velocityField->SetSpacing (m_Spacing);
    velocityField->SetOrigin (m_Origin);
    velocityField->SetDirection (m_Direction);
    velocityField->Allocate();

    itk::ImageRegionIteratorWithIndex < VelocityFieldType > svfIterator(velocityField,m_LargestRegion);
    tmpMatrix.fill(0);
    VelocityFieldIndexType curIndex;
    VelocityFieldPixelType curDisp;

    rigidIterator = itk::ImageRegionIterator < RigidFieldType > (rigidField,m_LargestRegion);
    while (!svfIterator.IsAtEnd())
    {
        curLog = rigidIterator.Get();

        unsigned int pos = 0;
        for (unsigned int j = 0;j < NDimensions;++j)
            for (unsigned int k = j+1;k < NDimensions;++k)
            {
                tmpMatrix(j,k) = curLog[pos];
                tmpMatrix(k,j) = - curLog[pos];
                ++pos;
            }

        for (unsigned int j = 0;j < NDimensions;++j)
        {
            tmpMatrix(j,NDimensions) = curLog[pos];
            ++pos;
        }

        curIndex = svfIterator.GetIndex();
        velocityField->TransformIndexToPhysicalPoint(curIndex,curPoint);
        for (unsigned int i = 0;i < NDimensions;++i)
        {
            curDisp[i] = tmpMatrix(i,NDimensions);
            for (unsigned int j = 0;j < NDimensions;++j)
                curDisp[i] += tmpMatrix(i,j) * curPoint[j];
        }

        svfIterator.Set(curDisp);

        ++rigidIterator;
        ++svfIterator;
    }

    // Create the final transform
    typename BaseOutputTransformType::Pointer resultTransform = BaseOutputTransformType::New();
    resultTransform->SetIdentity();
    resultTransform->SetParametersAsVectorField(velocityField.GetPointer());

    this->SetOutput(resultTransform);
}


template <unsigned int NDimensions>
void
BalooSVFTransformAgregator <NDimensions>::
estimateSVFFromAffineTransforms()
{
    const unsigned int NDegreesFreedom = NDimensions * (NDimensions + 1);
    typedef itk::Image < itk::Vector <ScalarType, NDegreesFreedom>, NDimensions > AffineFieldType;
    typedef typename AffineFieldType::Pointer AffineFieldPointer;
    typedef itk::Vector <ScalarType, NDegreesFreedom> AffineVectorType;

    unsigned int nbPts = this->GetInputRegions().size();
    AffineFieldPointer affineField = AffineFieldType::New();
    affineField->Initialize();

    affineField->SetRegions (m_LargestRegion);
    affineField->SetSpacing (m_Spacing);
    affineField->SetOrigin (m_Origin);
    affineField->SetDirection (m_Direction);
    affineField->Allocate();

    AffineVectorType zeroDisp;
    zeroDisp.Fill(0);
    itk::ImageRegionIterator < AffineFieldType > affineIterator(affineField,m_LargestRegion);
    while (!affineIterator.IsAtEnd())
    {
        affineIterator.Set(zeroDisp);
        ++affineIterator;
    }

    WeightImagePointer weights = WeightImageType::New();
    weights->Initialize();

    weights->SetRegions (m_LargestRegion);
    weights->SetSpacing (m_Spacing);
    weights->SetOrigin (m_Origin);
    weights->SetDirection (m_Direction);
    weights->Allocate();
    weights->FillBuffer(anima::zeroWeight);

    AffineVectorType curLog;
    VelocityFieldPointType curPoint;

    vnl_matrix <InternalScalarType> tmpMatrix(NDimensions+1,NDimensions+1,0);
    tmpMatrix(NDimensions,NDimensions) = 1;
    std::vector < AffineVectorType > logVectors(nbPts);
    std::vector < VelocityFieldIndexType > posIndexes(nbPts);

    if (dynamic_cast <anima::DirectionScaleSkewTransform <InternalScalarType> *> (this->GetInputTransform(0)) == 0)
    {
        typedef anima::MatrixLoggerFilter<InternalScalarType,ScalarType,NDimensions,NDegreesFreedom> MatrixLoggerFilterType;

        MatrixLoggerFilterType *logFilter = new MatrixLoggerFilterType;
        logFilter->SetInput(this->GetInputTransforms());
        logFilter->SetNumberOfThreads(m_NumberOfThreads);
        logFilter->SetUseRigidTransforms(false);

        logFilter->Update();
        logVectors = logFilter->GetOutput();

        delete logFilter;
    }
    else
    {
        typedef anima::DirectionScaleSkewTransform <InternalScalarType> DistoTrsfType;
        DistoTrsfType *tmpTrsf;
        for (unsigned int i = 0;i < nbPts;++i)
        {
            tmpTrsf = dynamic_cast <DistoTrsfType *> (this->GetInputTransform(i));
            logVectors[i] = tmpTrsf->GetLogVector();
        }
    }

    for (unsigned int i = 0;i < nbPts;++i)
    {
        double tmpWeight = this->GetInputWeight(i);
        weights->TransformPhysicalPointToIndex(this->GetInputOrigin(i),posIndexes[i]);

        if (boost::math::isnan(logVectors[i][0]))
        {
            logVectors[i].Fill(0);
            tmpWeight = 0;
        }

        affineField->SetPixel(posIndexes[i],logVectors[i] * tmpWeight);
        weights->SetPixel(posIndexes[i],tmpWeight);
    }

    for (unsigned int i = 0;i < m_DamIndexes.size();++i)
        weights->SetPixel(m_DamIndexes[i],1.0);

    anima::filterInputs<ScalarType,NDimensions,NDegreesFreedom>(weights,affineField,logVectors,posIndexes,this);

    VelocityFieldPointer velocityField = VelocityFieldType::New();
    velocityField->Initialize();

    velocityField->SetRegions (m_LargestRegion);
    velocityField->SetSpacing (m_Spacing);
    velocityField->SetOrigin (m_Origin);
    velocityField->SetDirection (m_Direction);
    velocityField->Allocate();

    itk::ImageRegionIteratorWithIndex < VelocityFieldType > svfIterator(velocityField,m_LargestRegion);
    tmpMatrix.fill(0);
    VelocityFieldIndexType curIndex;
    VelocityFieldPixelType curDisp;

    affineIterator = itk::ImageRegionIterator < AffineFieldType > (affineField,m_LargestRegion);
    while (!svfIterator.IsAtEnd())
    {
        curLog = affineIterator.Get();

        unsigned int pos = 0;
        for (unsigned int j = 0;j < NDimensions;++j)
            for (unsigned int k = 0;k <= NDimensions;++k)
            {
                tmpMatrix(j,k) = curLog[pos];
                ++pos;
            }

        curIndex = svfIterator.GetIndex();
        velocityField->TransformIndexToPhysicalPoint(curIndex,curPoint);
        for (unsigned int i = 0;i < NDimensions;++i)
        {
            curDisp[i] = tmpMatrix(i,NDimensions);
            for (unsigned int j = 0;j < NDimensions;++j)
                curDisp[i] += tmpMatrix(i,j) * curPoint[j];
        }

        svfIterator.Set(curDisp);

        ++affineIterator;
        ++svfIterator;
    }

    // Create the final transform
    typename BaseOutputTransformType::Pointer resultTransform = BaseOutputTransformType::New();
    resultTransform->SetIdentity();
    resultTransform->SetParametersAsVectorField(velocityField.GetPointer());

    this->SetOutput(resultTransform);
}

//---------------------------------------------------------------------------------------------
// global filterInputs funtion
//---------------------------------------------------------------------------------------------

template <class ScalarType, unsigned int NDimensions, unsigned int NDegreesOfFreedom>
void
filterInputs(itk::Image <ScalarType,NDimensions> *weights,
             typename itk::Image < itk::Vector <ScalarType, NDegreesOfFreedom>, NDimensions >::Pointer &output,
             std::vector < itk::Vector <ScalarType, NDegreesOfFreedom> > &curTrsfs,
             std::vector < typename itk::Image < itk::Vector <ScalarType, NDegreesOfFreedom>, NDimensions >::IndexType > &posIndexes,
             BalooSVFTransformAgregator <NDimensions> *filterPtr)
{
    typedef itk::Image < itk::Vector <ScalarType, NDegreesOfFreedom>, NDimensions > FieldType;
    typedef itk::Vector <ScalarType, NDegreesOfFreedom> FieldPixelType;

    typedef itk::Image <ScalarType,NDimensions> WeightImageType;

    FieldPixelType zeroTrsf;
    zeroTrsf.Fill(0);

    typedef anima::SmoothingRecursiveYvvGaussianImageFilter<WeightImageType,WeightImageType> WeightSmootherType;
    typename WeightSmootherType::Pointer weightSmooth = WeightSmootherType::New();

    weightSmooth->SetInput(weights);
    weightSmooth->SetSigma(filterPtr->GetExtrapolationSigma());
    weightSmooth->SetNumberOfThreads(filterPtr->GetNumberOfThreads());

    weightSmooth->Update();

    typename WeightImageType::Pointer smoothedWeights = weightSmooth->GetOutput();
    smoothedWeights->DisconnectPipeline();
    weightSmooth = 0;

    typedef anima::SmoothingRecursiveYvvGaussianImageFilter <FieldType, FieldType> SVFSmoothingFilterType;
    typename SVFSmoothingFilterType::Pointer smootherPtr = SVFSmoothingFilterType::New();

    smootherPtr->SetInput(output);
    smootherPtr->SetSigma(filterPtr->GetExtrapolationSigma());
    smootherPtr->SetNumberOfThreads(filterPtr->GetNumberOfThreads());

    smootherPtr->Update();

    typename FieldType::Pointer smoothedField = smootherPtr->GetOutput();
    smoothedField->DisconnectPipeline();
    smootherPtr = 0;

    // Now remove outliers
    FieldPixelType curTrsf;

    unsigned int nbPts = filterPtr->GetInputRegions().size();
    std::vector < double > residuals(nbPts);
    double averageResidual = 0, varResidual = 0;
    for (unsigned int i = 0;i < nbPts;++i)
    {
        double weightSmoothed = smoothedWeights->GetPixel(posIndexes[i]);
        if (weightSmoothed > 0)
        {
            curTrsf = smoothedField->GetPixel(posIndexes[i]);
            curTrsf /= weightSmoothed;
        }
        else
            curTrsf = zeroTrsf;

        double residual = 0;

        for (unsigned int j = 0;j < NDegreesOfFreedom;++j)
            residual += (curTrsf[j] - curTrsfs[i][j]) * (curTrsf[j] - curTrsfs[i][j]);

        varResidual += residual;
        residual = sqrt(residual);
        averageResidual += residual;
        residuals[i] = residual;
    }

    smoothedField = 0;
    smoothedWeights = 0;

    varResidual = sqrt((varResidual - averageResidual * averageResidual / nbPts) / (nbPts - 1.0));
    averageResidual /= nbPts;

    for (unsigned int i = 0;i < nbPts;++i)
    {
        if (residuals[i] > averageResidual + filterPtr->GetOutlierRejectionSigma() * varResidual)
        {
            output->SetPixel(posIndexes[i],zeroTrsf);
            weights->SetPixel(posIndexes[i],anima::zeroWeight);
        }
    }

    weightSmooth = WeightSmootherType::New();

    weightSmooth->SetInput(weights);
    weightSmooth->SetSigma(filterPtr->GetExtrapolationSigma());
    weightSmooth->SetNumberOfThreads(filterPtr->GetNumberOfThreads());

    weightSmooth->Update();

    smoothedWeights = weightSmooth->GetOutput();
    smoothedWeights->DisconnectPipeline();

    smootherPtr = SVFSmoothingFilterType::New();

    smootherPtr->SetInput(output);
    smootherPtr->SetSigma(filterPtr->GetExtrapolationSigma());
    smootherPtr->SetNumberOfThreads(filterPtr->GetNumberOfThreads());

    smootherPtr->Update();

    output = smootherPtr->GetOutput();
    output->DisconnectPipeline();

    typename WeightImageType::RegionType largestRegion = weights->GetLargestPossibleRegion();
    itk::ImageRegionIterator < FieldType > svfIterator = itk::ImageRegionIterator < FieldType > (output,largestRegion);
    itk::ImageRegionConstIterator < WeightImageType > weightItr(smoothedWeights,largestRegion);

    while (!svfIterator.IsAtEnd())
    {
        if (weightItr.Value() > 0)
        {
            curTrsf = svfIterator.Get();
            curTrsf /= weightItr.Value();
            svfIterator.Set(curTrsf);
        }
        else
            svfIterator.Set(zeroTrsf);

        ++weightItr;
        ++svfIterator;
    }
}

} // end of namespace anima
