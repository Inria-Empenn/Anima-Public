#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkInterpolateImageFunction.h>

#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template <typename TImageType, typename TInterpolatorPrecisionType=float>
class OrientedModelBaseResampleImageFilter :
        public itk::ImageToImageFilter <TImageType, TImageType>
{
public:
    /** Standard class typedefs. */
    typedef OrientedModelBaseResampleImageFilter Self;
    typedef TImageType InputImageType;
    typedef TImageType TOutputImage;
    itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

    typedef itk::Image <unsigned char, itkGetStaticConstMacro(ImageDimension)> GeometryImageType;
    typedef typename GeometryImageType::Pointer GeometryImagePointer;

    typedef itk::ImageToImageFilter<TImageType, TImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef itk::InterpolateImageFunction<InputImageType, TInterpolatorPrecisionType> InterpolatorType;
    typedef typename InterpolatorType::Pointer  InterpolatorPointer;

    typedef typename InterpolatorType::ContinuousIndexType ContinuousIndexType;
    typedef typename InterpolatorType::PointType PointType;

    /** Run-time type information (and related methods) */
    itkTypeMacro(OrientedModelBaseResampleImageFilter, ImageToImageFilter);

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename InputImageType::IndexType InputIndexType;

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename TOutputImage::Pointer OutputImagePointer;
    typedef typename TOutputImage::RegionType OutputImageRegionType;

    typedef itk::Transform <TInterpolatorPrecisionType,
    itkGetStaticConstMacro(ImageDimension),
    itkGetStaticConstMacro(ImageDimension)> TransformType;

    typedef itk::MatrixOffsetTransformBase <TInterpolatorPrecisionType,
    itkGetStaticConstMacro(ImageDimension),
    itkGetStaticConstMacro(ImageDimension)> MatrixTransformType;

    typedef typename TransformType::Pointer TransformPointer;

    void SetTransform (TransformType *trsf)
    {
        m_Transform = trsf;
        if (dynamic_cast <MatrixTransformType *> (m_Transform.GetPointer()) != 0)
            m_LinearTransform = true;
        else
            m_LinearTransform = false;
    }

    itkGetMacro(Transform,TransformType *);

    itkSetMacro(Interpolator,InterpolatorPointer);
    itkGetMacro(Interpolator,InterpolatorType *);

    typedef typename TOutputImage::SpacingType   SpacingType;
    typedef typename TOutputImage::PointType     OriginPointType;
    typedef typename TOutputImage::DirectionType DirectionType;
    typedef typename TOutputImage::RegionType RegionType;

    itkSetMacro( OutputSpacing, SpacingType );
    itkGetConstReferenceMacro( OutputSpacing, SpacingType );

    itkSetMacro( OutputOrigin, OriginPointType );
    itkGetConstReferenceMacro( OutputOrigin, OriginPointType );

    itkSetMacro( OutputDirection, DirectionType );
    itkGetConstReferenceMacro( OutputDirection, DirectionType );

    itkSetMacro( OutputLargestPossibleRegion, RegionType );
    itkGetConstReferenceMacro( OutputLargestPossibleRegion, RegionType );

protected:
    OrientedModelBaseResampleImageFilter()
    {
        m_LinearTransform = false;

        m_Transform = 0;
        m_Interpolator = 0;
    }

    virtual ~OrientedModelBaseResampleImageFilter() {}

    virtual void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);

    virtual unsigned int GetOutputVectorLength();

    bool isZero(const InputPixelType &dataVec)
    {
        unsigned int vectorSize = dataVec.GetNumberOfElements();

        for (unsigned int i = 0;i < vectorSize;++i)
        {
            if (dataVec[i] != 0)
                return false;
        }

        return true;
    }

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

    void LinearThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);
    void NonLinearThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);

    vnl_matrix <double> ComputeLinearRotationMatrix();
    void ComputeLocalRotationMatrix(InputIndexType &index, vnl_matrix <double> &rotMatrix);

    //! Initializes the default interpolator, might change in derived classes
    virtual void InitializeInterpolator();

    //! For some models, the transform induced rotation matrix is not enough and needs to be used to get rotation parameters, default is just a copy of the current matrix
    virtual void ComputeRotationParametersFromRotationMatrix(const vnl_matrix <double> &transformRotationMatrix, vnl_matrix <double> &modelRotationMatrix)
    {
        modelRotationMatrix = transformRotationMatrix;
    }

    //! Needs to be implemented in sub-classes, does the actual rotation of the model
    virtual void RotateInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelRotationMatrix,
                                         InputPixelType &rotatedModel, itk::ThreadIdType threadId) = 0;

private:
    OrientedModelBaseResampleImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    TransformPointer m_Transform;
    bool m_LinearTransform;

    InterpolatorPointer m_Interpolator;

    SpacingType             m_OutputSpacing;
    OriginPointType         m_OutputOrigin;
    DirectionType           m_OutputDirection;
    RegionType              m_OutputLargestPossibleRegion;

    InputIndexType m_StartIndex, m_EndIndex;
    InputIndexType m_StartIndDef, m_EndIndDef;
};

} // end namespace anima

#include "animaOrientedModelBaseResampleImageFilter.hxx"
