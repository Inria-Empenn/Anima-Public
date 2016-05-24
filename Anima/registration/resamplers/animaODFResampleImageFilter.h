#pragma once

#include <animaOrientedModelBaseResampleImageFilter.h>

#include <vector>
#include <vnl/vnl_matrix.h>

namespace anima
{
template <typename TImageType, typename TInterpolatorPrecisionType=float>
class ODFResampleImageFilter :
        public OrientedModelBaseResampleImageFilter <TImageType,TInterpolatorPrecisionType>
{
public:
    /** Standard class typedefs. */
    typedef ODFResampleImageFilter Self;

    typedef OrientedModelBaseResampleImageFilter <TImageType,TInterpolatorPrecisionType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef typename Superclass::InputPixelType InputPixelType;
    typedef typename Superclass::InputImageType InputImageType;
    itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(ODFResampleImageFilter, OrientedModelBaseResampleImageFilter)

protected:
    ODFResampleImageFilter()
    {
        m_LOrder = 4;
    }

    virtual ~ODFResampleImageFilter() {}

    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

    virtual void ReorientInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelOrientationMatrix,
                                           InputPixelType &rotatedModel, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ODFResampleImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_LOrder;

    std::vector < std::vector <double> > m_EulerAngles;
    std::vector < vnl_matrix <double> > m_ODFRotationMatrices;
};

} // end namespace anima

#include "animaODFResampleImageFilter.hxx"
