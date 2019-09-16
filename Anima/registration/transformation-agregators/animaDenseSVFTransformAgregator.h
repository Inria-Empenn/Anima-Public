#pragma once

#include <animaBaseTransformAgregator.h>
#include <itkStationaryVelocityFieldTransform.h>

namespace anima
{

template <unsigned int NDimensions = 3>
class DenseSVFTransformAgregator :
public BaseTransformAgregator <NDimensions>
{
public:
    typedef BaseTransformAgregator <NDimensions> Superclass;
    typedef typename Superclass::PointType PointType;
    typedef typename Superclass::BaseInputTransformType BaseInputTransformType;
    typedef typename Superclass::ScalarType ScalarType;
    typedef typename itk::StationaryVelocityFieldTransform<ScalarType,NDimensions> BaseOutputTransformType;
    typedef typename Superclass::InternalScalarType InternalScalarType;

    typedef itk::Image <ScalarType,NDimensions> WeightImageType;
    typedef typename WeightImageType::IndexType IndexType;
    typedef typename WeightImageType::Pointer WeightImagePointer;

    typedef itk::MatrixOffsetTransformBase <InternalScalarType, NDimensions, NDimensions> BaseMatrixTransformType;
    typedef typename BaseMatrixTransformType::ParametersType ParametersType;

    typedef typename BaseOutputTransformType::VectorFieldType VelocityFieldType;
    typedef typename VelocityFieldType::Pointer VelocityFieldPointer;
    typedef typename VelocityFieldType::IndexType VelocityFieldIndexType;
    typedef typename VelocityFieldType::PixelType VelocityFieldPixelType;
    typedef typename VelocityFieldType::PointType VelocityFieldPointType;
    typedef typename VelocityFieldType::RegionType VelocityFieldRegionType;
    typedef typename VelocityFieldType::SpacingType VelocityFieldSpacingType;
    typedef typename VelocityFieldType::DirectionType VelocityFieldDirectionType;

    DenseSVFTransformAgregator();
    virtual ~DenseSVFTransformAgregator() {}

    virtual bool Update();

    void SetExtrapolationSigma(double sigma) {m_ExtrapolationSigma = sigma;}
    void SetOutlierRejectionSigma(double sigma) {m_OutlierRejectionSigma = sigma;}

    void SetNumberOfWorkUnits(unsigned int num) {m_NumberOfThreads = num;}

    void SetNeighborhoodHalfSize(unsigned int num) {m_NeighborhoodHalfSize = num;}
    void SetDistanceBoundary(double num) {m_DistanceBoundary = num;}
    void SetMEstimateConvergenceThreshold(double num) {m_MEstimateConvergenceThreshold = num;}

    template <class TInputImageType> void SetGeometryInformation(const TInputImageType *geomImage)
    {
        if (geomImage == NULL)
            return;

        m_Origin = geomImage->GetOrigin();
        m_LargestRegion = geomImage->GetLargestPossibleRegion();
        m_Spacing = geomImage->GetSpacing();
        m_Direction = geomImage->GetDirection();
    }

protected:
    double m_ExtrapolationSigma;
    double m_OutlierRejectionSigma;

    unsigned int m_NeighborhoodHalfSize;
    double m_DistanceBoundary;
    double m_MEstimateConvergenceThreshold;

    VelocityFieldPointType m_Origin;
    VelocityFieldRegionType m_LargestRegion;
    VelocityFieldSpacingType m_Spacing;
    VelocityFieldDirectionType m_Direction;

    unsigned int m_NumberOfThreads;

private:
    void estimateSVFFromTranslations();
    void estimateSVFFromRigidTransforms();
    void estimateSVFFromAffineTransforms();
};

} // end of namespace anima

#include "animaDenseSVFTransformAgregator.hxx"
