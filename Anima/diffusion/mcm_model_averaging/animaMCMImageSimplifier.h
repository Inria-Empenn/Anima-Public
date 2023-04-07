#pragma once
#include <animaNumberedThreadImageToImageFilter.h>
#include <animaMCMImage.h>
#include <itkImage.h>
#include <itkSymmetricEigenAnalysis.h>

#include <animaBaseCompartment.h>
#include <animaMultiCompartmentModel.h>

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

namespace anima
{

template <class PixelScalarType>
class MCMImageSimplifier :
public anima::NumberedThreadImageToImageFilter < anima::MCMImage<PixelScalarType,3>, anima::MCMImage<PixelScalarType,3> >
{
public:
    /** Standard class typedefs. */
    typedef MCMImageSimplifier Self;
    typedef anima::MCMImage<PixelScalarType,3> InputImageType;
    typedef anima::MCMImage<PixelScalarType,3> OutputImageType;
    typedef itk::Image<unsigned int,3> MoseImageType;
    typedef anima::NumberedThreadImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef typename MCModelType::Pointer MCModelPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMImageSimplifier, anima::NumberedThreadImageToImageFilter)

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename OutputImageType::PixelType OutputPixelType;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetMoseVolume(MoseImageType *vol) {m_MoseMap = vol;}

protected:
    MCMImageSimplifier()
    : Superclass()
    {
        m_MoseMap = 0;
    }

    virtual ~MCMImageSimplifier() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    void InitializeReferenceOutputModel();

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMImageSimplifier);

    MoseImageType::Pointer m_MoseMap;
    MCModelPointer m_ReferenceOutputModel;
};
    
} // end of namespace anima

#include "animaMCMImageSimplifier.hxx"
