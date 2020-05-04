#pragma once

#include <animaOrientedModelBaseResampleImageFilter.h>
#include <animaMultiCompartmentModel.h>

namespace anima
{

template <typename TImageType, typename TInterpolatorPrecisionType=double>
class MCMResampleImageFilter :
        public OrientedModelBaseResampleImageFilter <TImageType, TInterpolatorPrecisionType>
{
public:
    /** Standard class typedefs. */
    typedef MCMResampleImageFilter Self;

    typedef OrientedModelBaseResampleImageFilter <TImageType,TInterpolatorPrecisionType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef typename Superclass::InputPixelType InputPixelType;
    typedef typename Superclass::InputImageType InputImageType;
    itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMResampleImageFilter, OrientedModelBaseResampleImageFilter)

    //! Sets reference output MCM model, necessary to determine output organization (and rotate)
    void SetReferenceOutputModel(const MCModelPointer &model);

protected:
    MCMResampleImageFilter()
    {
    }

    virtual ~MCMResampleImageFilter() {}

    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
    virtual void InitializeInterpolator() ITK_OVERRIDE;

    virtual void ReorientInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelOrientationMatrix,
                                           InputPixelType &orientedModel, itk::ThreadIdType threadId) ITK_OVERRIDE;

    virtual unsigned int GetOutputVectorLength() ITK_OVERRIDE;

    virtual itk::LightObject::Pointer InternalClone() const ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMResampleImageFilter);

    MCModelPointer m_ReferenceOutputModel;

    std::vector <MCModelPointer> m_WorkModels;
};

} // end namespace anima

#include "animaMCMResampleImageFilter.hxx"
