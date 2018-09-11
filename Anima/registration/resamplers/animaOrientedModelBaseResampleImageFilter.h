#pragma once

#include <iostream>
#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkInterpolateImageFunction.h>

#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template <typename TImageType, typename TInterpolatorPrecisionType=float>
class OrientedModelBaseResampleImageFilter :
        public anima::MaskedImageToImageFilter <TImageType, TImageType>
{
public:
    /** Standard class typedefs. */
    typedef OrientedModelBaseResampleImageFilter Self;
    typedef TImageType InputImageType;
    typedef TImageType TOutputImage;
    itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

    typedef itk::Image <unsigned char, itkGetStaticConstMacro(ImageDimension)> GeometryImageType;
    typedef typename GeometryImageType::Pointer GeometryImagePointer;

    typedef anima::MaskedImageToImageFilter<TImageType, TImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef itk::InterpolateImageFunction<InputImageType, TInterpolatorPrecisionType> InterpolatorType;
    typedef typename InterpolatorType::Pointer  InterpolatorPointer;

    typedef typename InterpolatorType::ContinuousIndexType ContinuousIndexType;
    typedef typename InterpolatorType::PointType PointType;

    /** Run-time type information (and related methods) */
    itkTypeMacro(OrientedModelBaseResampleImageFilter, anima::MaskedImageToImageFilter)

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename InputImageType::PointType InputPointType;
    typedef typename InputImageType::IndexType InputIndexType;

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename TOutputImage::Pointer OutputImagePointer;
    typedef typename TOutputImage::RegionType OutputImageRegionType;

    typedef itk::Transform <TInterpolatorPrecisionType,
    itkGetStaticConstMacro(ImageDimension),
    itkGetStaticConstMacro(ImageDimension)> TransformType;

    typedef typename TransformType::Pointer TransformPointer;

    itkSetObjectMacro(Transform,TransformType)
    itkGetObjectMacro(Transform,TransformType)

    itkSetMacro(Interpolator,InterpolatorPointer)
    itkGetMacro(Interpolator,InterpolatorType *)

    itkGetMacro(FiniteStrainReorientation, bool)
    itkSetMacro(FiniteStrainReorientation, bool)

    typedef typename TOutputImage::SpacingType SpacingType;
    typedef typename TOutputImage::PointType OriginPointType;
    typedef typename TOutputImage::DirectionType DirectionType;
    typedef typename TOutputImage::RegionType RegionType;

    itkSetMacro(OutputSpacing, SpacingType)
    itkGetConstReferenceMacro(OutputSpacing, SpacingType)

    itkSetMacro(OutputOrigin, OriginPointType)
    itkGetConstReferenceMacro(OutputOrigin, OriginPointType)

    itkSetMacro(OutputDirection, DirectionType)
    itkGetConstReferenceMacro(OutputDirection, DirectionType)

    itkSetMacro(OutputLargestPossibleRegion, RegionType)
    itkGetConstReferenceMacro(OutputLargestPossibleRegion, RegionType)

protected:
    OrientedModelBaseResampleImageFilter()
    {
        m_Transform = 0;
        m_Interpolator = 0;

        m_FiniteStrainReorientation = true;
    }

    virtual ~OrientedModelBaseResampleImageFilter() {}

    virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;
    virtual void GenerateOutputInformation() ITK_OVERRIDE;

    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    virtual unsigned int GetOutputVectorLength();

    template <class T>
    inline bool isZero(itk::VariableLengthVector <T> &dataVec)
    {
        unsigned int vectorSize = dataVec.GetNumberOfElements();

        for (unsigned int i = 0;i < vectorSize;++i)
        {
            if (dataVec[i] != 0)
                return false;
        }

        return true;
    }

    // Fake method for compilation purposes, should never go in there
    template <class T> inline bool isZero(T &data)
    {
        itkExceptionMacro("Access to unauthorized method");
        return true;
    }

    //! Utility function to initialize output images pixel to zero for vector images
    template <class T> inline void InitializeZeroPixel(itk::VariableLengthVector <T> &zeroPixel)
    {
        zeroPixel.Fill(0.0);
    }

    //! Utility function to initialize output images pixel to zero for all images except vector images
    template <class T> inline void InitializeZeroPixel(T &zeroPixel)
    {
        zeroPixel = itk::NumericTraits <T>::ZeroValue();
    }

    void LinearThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, unsigned int threadId);
    void NonLinearThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, unsigned int threadId);

    /**
     * Compute local re-orientation transformation matrix from linear transform
     * (rotation if finite strain, Jacobian otherwise)
     */
    vnl_matrix <double> ComputeLinearJacobianMatrix();

    /**
     * Compute local re-orientation transformation matrix from non linear transform
     * (rotation if finite strain, Jacobian otherwise)
     */
    void ComputeLocalJacobianMatrix(InputIndexType &index, vnl_matrix <double> &reorientationMatrix);

    //! Initializes the default interpolator, might change in derived classes
    virtual void InitializeInterpolator();

    /**
     * For some models / re-orientation schemes, the transform induced Jacobiam matrix needs to be used to get rotation parameters
     * if needed, this method should be overloaded
     */
    virtual void ComputeRotationParametersFromReorientationMatrix(vnl_matrix <double> &reorientationMatrix,
                                                                  vnl_matrix <double> &modelOrientationMatrix);

    //! Needs to be implemented in sub-classes, does the actual re-orientation of the model
    virtual void ReorientInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelOrientationMatrix,
                                           InputPixelType &orientedModel, itk::ThreadIdType threadId) = 0;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(OrientedModelBaseResampleImageFilter);

    TransformPointer m_Transform;
    bool m_FiniteStrainReorientation;
    InterpolatorPointer m_Interpolator;

    SpacingType m_OutputSpacing;
    OriginPointType m_OutputOrigin;
    DirectionType m_OutputDirection;
    RegionType m_OutputLargestPossibleRegion;

    InputIndexType m_StartIndDef, m_EndIndDef;
};

} // end namespace anima

#include "animaOrientedModelBaseResampleImageFilter.hxx"
