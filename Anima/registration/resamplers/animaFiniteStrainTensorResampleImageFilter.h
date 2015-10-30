#pragma once

#include <animaOrientedModelBaseResampleImageFilter.h>

namespace anima
{

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType=float>
class FiniteStrainTensorResampleImageFilter :
        public OrientedModelBaseResampleImageFilter <TInputScalarType,Dimension,TInterpolatorPrecisionType>
{
public:
    /** Standard class typedefs. */
    typedef FiniteStrainTensorResampleImageFilter Self;

    typedef OrientedModelBaseResampleImageFilter <TInputScalarType,Dimension,TInterpolatorPrecisionType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef typename Superclass::InputPixelType InputPixelType;
    typedef typename Superclass::InputImageType InputImageType;
    itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(FiniteStrainTensorResampleImageFilter, OrientedModelBaseResampleImageFilter);

protected:
    FiniteStrainTensorResampleImageFilter()
    {
        m_VectorSize = ImageDimension * (ImageDimension + 1) / 2;
        m_TensorDimension = ImageDimension;
    }

    virtual ~FiniteStrainTensorResampleImageFilter() {}

    virtual void RotateInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelRotationMatrix,
                                         InputPixelType &rotatedModel);

private:
    FiniteStrainTensorResampleImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_VectorSize;
    unsigned int m_TensorDimension;
};

} // end namespace anima

#include "animaFiniteStrainTensorResampleImageFilter.hxx"
