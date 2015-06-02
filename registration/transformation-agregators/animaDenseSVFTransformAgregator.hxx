#pragma once

#include <animaMEstimateSVFImageFilter.h>
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
DenseSVFTransformAgregator <NDimensions>::
DenseSVFTransformAgregator() : Superclass()
{
    m_ExtrapolationSigma = 4.0;
    m_OutlierRejectionSigma = 3.0;

    m_NeighborhoodHalfSize = (unsigned int)floor(m_ExtrapolationSigma * 3);
    m_DistanceBoundary = m_ExtrapolationSigma * 3;
    m_MEstimateConvergenceThreshold = 0.001;

    m_NumberOfThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
}

template <unsigned int NDimensions>
bool
DenseSVFTransformAgregator <NDimensions>::
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
DenseSVFTransformAgregator <NDimensions>::
estimateSVFFromTranslations()
{
    const unsigned int NDegreesFreedom = NDimensions;
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
    weights->FillBuffer(0);

    ParametersType tmpParams;
    VelocityFieldPixelType curDisp;
    VelocityFieldIndexType posIndex;
    VelocityFieldRegionType tmpRegion;

    for (unsigned int i = 0;i < nbPts;++i)
    {
        tmpParams = this->GetInputTransform(i)->GetParameters();
        tmpRegion = this->GetInputRegions()[i];
        double tmpWeight = this->GetInputWeight(i);

        for (unsigned int j = 0;j < NDimensions;++j)
        {
            //round() is not always defined
            //round() can be replaced by floor(i+0.5)
            posIndex[j] = (unsigned int)floor((tmpRegion.GetIndex()[j] + tmpRegion.GetSize()[j] / 2.0) + 0.5);
            curDisp[j] = tmpParams[j];
        }

        velocityField->SetPixel(posIndex,curDisp);
        weights->SetPixel(posIndex,tmpWeight);
    }

    for (unsigned int i = 0;i < m_DamIndexes.size();++i)
        weights->SetPixel(m_DamIndexes[i],1.0);

    typedef anima::MEstimateSVFImageFilter <ScalarType,NDegreesFreedom,NDimensions> SVFMEstimateType;
    typename SVFMEstimateType::Pointer fieldSmoother = SVFMEstimateType::New();

    fieldSmoother->SetInput(velocityField);
    fieldSmoother->SetWeightImage(weights);
    fieldSmoother->SetFluidSigma(m_ExtrapolationSigma);
    fieldSmoother->SetNeighborhoodHalfSize(m_NeighborhoodHalfSize);
    fieldSmoother->SetDistanceBoundary(m_DistanceBoundary);
    fieldSmoother->SetMEstimateFactor(m_OutlierRejectionSigma);

    fieldSmoother->SetConvergenceThreshold(m_MEstimateConvergenceThreshold);
    fieldSmoother->SetMaxNumIterations(100);

    fieldSmoother->SetNumberOfThreads(m_NumberOfThreads);

    fieldSmoother->Update();

    velocityField = fieldSmoother->GetOutput();
    velocityField->DisconnectPipeline();

    // Create the final transform
    typename BaseOutputTransformType::Pointer resultTransform = BaseOutputTransformType::New();
    resultTransform->SetIdentity();
    resultTransform->SetParametersAsVectorField(velocityField.GetPointer());

    this->SetOutput(resultTransform);
}


template <unsigned int NDimensions>
void
DenseSVFTransformAgregator <NDimensions>::
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
    weights->FillBuffer(0);

    RigidVectorType curLog;
    VelocityFieldIndexType posIndex;
    VelocityFieldRegionType tmpRegion;

    vnl_matrix <InternalScalarType> tmpMatrix(NDimensions+1,NDimensions+1,0);
    tmpMatrix(NDimensions,NDimensions) = 1;

    typedef anima::LogRigid3DTransform <InternalScalarType> LogRigidTransformType;

    for (unsigned int i = 0;i < nbPts;++i)
    {
        tmpRegion = this->GetInputRegions()[i];
        double tmpWeight = this->GetInputWeight(i);

        for (unsigned int j = 0;j < NDimensions;++j)
            posIndex[j] = (unsigned int)floor((tmpRegion.GetIndex()[j] + tmpRegion.GetSize()[j] / 2.0) + 0.5);

        LogRigidTransformType *tmpTrsf = (LogRigidTransformType *)this->GetInputTransform(i);
        rigidField->SetPixel(posIndex,tmpTrsf->GetLogVector());
        weights->SetPixel(posIndex,tmpWeight);
    }

    for (unsigned int i = 0;i < m_DamIndexes.size();++i)
        weights->SetPixel(m_DamIndexes[i],1.0);

    typedef anima::MEstimateSVFImageFilter <ScalarType,NDegreesFreedom,NDimensions> SVFMEstimateType;
    typename SVFMEstimateType::Pointer fieldSmoother = SVFMEstimateType::New();

    fieldSmoother->SetInput(rigidField);
    fieldSmoother->SetWeightImage(weights);
    fieldSmoother->SetFluidSigma(m_ExtrapolationSigma);
    fieldSmoother->SetNeighborhoodHalfSize(m_NeighborhoodHalfSize);
    fieldSmoother->SetDistanceBoundary(m_DistanceBoundary);
    fieldSmoother->SetMEstimateFactor(m_OutlierRejectionSigma);

    fieldSmoother->SetConvergenceThreshold(m_MEstimateConvergenceThreshold);
    fieldSmoother->SetMaxNumIterations(100);

    fieldSmoother->SetNumberOfThreads(m_NumberOfThreads);

    fieldSmoother->Update();

    rigidField = fieldSmoother->GetOutput();
    rigidField->DisconnectPipeline();

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
    VelocityFieldPointType curPoint;
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
DenseSVFTransformAgregator <NDimensions>::
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
    weights->FillBuffer(0);

    AffineVectorType curLog;
    VelocityFieldIndexType posIndex;
    VelocityFieldRegionType tmpRegion;

    vnl_matrix <InternalScalarType> tmpMatrix(NDimensions+1,NDimensions+1,0);
    tmpMatrix(NDimensions,NDimensions) = 1;
    std::vector < AffineVectorType > logVectors(nbPts);

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
        tmpRegion = this->GetInputRegions()[i];
        double tmpWeight = this->GetInputWeight(i);

        for (unsigned int j = 0;j < NDimensions;++j)
            posIndex[j] = (unsigned int)floor((tmpRegion.GetIndex()[j] + tmpRegion.GetSize()[j] / 2.0) + 0.5);

        if (boost::math::isnan(logVectors[i][0]))
        {
            tmpWeight = 0;
            logVectors[i].Fill(0);
        }

        affineField->SetPixel(posIndex,logVectors[i]);
        weights->SetPixel(posIndex,tmpWeight);
    }

    for (unsigned int i = 0;i < m_DamIndexes.size();++i)
        weights->SetPixel(m_DamIndexes[i],1.0);

    typedef anima::MEstimateSVFImageFilter <ScalarType,NDegreesFreedom,NDimensions> SVFMEstimateType;
    typename SVFMEstimateType::Pointer fieldSmoother = SVFMEstimateType::New();

    fieldSmoother->SetInput(affineField);
    fieldSmoother->SetWeightImage(weights);
    fieldSmoother->SetFluidSigma(m_ExtrapolationSigma);
    fieldSmoother->SetNeighborhoodHalfSize(m_NeighborhoodHalfSize);
    fieldSmoother->SetDistanceBoundary(m_DistanceBoundary);
    fieldSmoother->SetMEstimateFactor(m_OutlierRejectionSigma);

    fieldSmoother->SetConvergenceThreshold(m_MEstimateConvergenceThreshold);
    fieldSmoother->SetMaxNumIterations(100);

    fieldSmoother->SetNumberOfThreads(m_NumberOfThreads);

    fieldSmoother->Update();

    affineField = fieldSmoother->GetOutput();
    affineField->DisconnectPipeline();

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
    VelocityFieldPointType curPoint;
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

} // end of namespace anima
