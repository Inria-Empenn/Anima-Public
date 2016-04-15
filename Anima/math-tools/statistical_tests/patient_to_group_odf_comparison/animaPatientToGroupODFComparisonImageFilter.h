#pragma once

#include <animaPatientToGroupComparisonImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

namespace anima
{

template <class PixelScalarType>
class PatientToGroupODFComparisonImageFilter :
        public PatientToGroupComparisonImageFilter <PixelScalarType>
{
public:
    /** Standard class typedefs. */
    typedef anima::PatientToGroupODFComparisonImageFilter<PixelScalarType> Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef anima::PatientToGroupComparisonImageFilter<PixelScalarType> Superclass;
    typedef typename Superclass::VectorType VectorType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(PatientToGroupODFComparisonImageFilter, PatientToGroupComparisonImageFilter);

    void SetSampleDirections (std::vector < std::vector < double > > &input) {m_SampleDirections = input;}

protected:
    PatientToGroupODFComparisonImageFilter()
        : Superclass()
    {
        m_SampleDirections.clear();
        m_ShData = NULL;
    }

    virtual ~PatientToGroupODFComparisonImageFilter()
    {
        if (m_ShData)
            delete m_ShData;
    }

    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
    unsigned int SampleFromDiffusionModels(std::vector <VectorType> &databaseValues, VectorType &patientVectorValue) ITK_OVERRIDE;

private:
    PatientToGroupODFComparisonImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector < std::vector <double> > m_SampleDirections;
    unsigned int m_LOrder;

    anima::ODFSphericalHarmonicBasis *m_ShData;
};

} // end namespace anima

#include "animaPatientToGroupODFComparisonImageFilter.hxx"
