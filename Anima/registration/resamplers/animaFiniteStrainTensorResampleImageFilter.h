#pragma once

#include <animaOrientedModelBaseResampleImageFilter.h>

namespace anima
{

template <typename TImageType, typename TInterpolatorPrecisionType=float>
class FiniteStrainTensorResampleImageFilter :
        public OrientedModelBaseResampleImageFilter <TImageType,TInterpolatorPrecisionType>
{
public:
    /** Standard class typedefs. */
    typedef FiniteStrainTensorResampleImageFilter Self;

    typedef OrientedModelBaseResampleImageFilter <TImageType,TInterpolatorPrecisionType> Superclass;
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

    virtual void BeforeThreadedGenerateData();
    virtual void RotateInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelRotationMatrix,
                                         InputPixelType &rotatedModel, itk::ThreadIdType threadId);

private:
    FiniteStrainTensorResampleImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_VectorSize;
    unsigned int m_TensorDimension;

    // Work variables
    std::vector < vnl_matrix <double> > m_WorkMats;
    std::vector < vnl_matrix <double> > m_TmpTensors;
};

} // end namespace anima

#include "animaFiniteStrainTensorResampleImageFilter.hxx"
