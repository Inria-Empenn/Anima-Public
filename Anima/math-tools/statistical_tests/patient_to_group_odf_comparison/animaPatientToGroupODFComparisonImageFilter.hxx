#pragma once
#include "animaPatientToGroupODFComparisonImageFilter.h"

namespace anima
{

template <class PixelScalarType>
void
PatientToGroupODFComparisonImageFilter<PixelScalarType>
::BeforeThreadedGenerateData ()
{
    Superclass::BeforeThreadedGenerateData();

    unsigned int ndim = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    m_LOrder = (unsigned int)floor((-3.0 + sqrt(8.0 * ndim + 1.0))/2.0);

    if (m_ShData)
        delete m_ShData;

    m_ShData = new anima::ODFSphericalHarmonicBasis (m_LOrder);
}

template <class PixelScalarType>
unsigned int
PatientToGroupODFComparisonImageFilter<PixelScalarType>
::SampleFromDiffusionModels(std::vector <VectorType> &databaseValues,
                            VectorType &patientVectorValue)
{
    if (m_SampleDirections.size() == 0)
        return patientVectorValue.GetSize();

    unsigned int numItems = databaseValues.size();
    for (unsigned int i = 0;i < numItems;++i)
        databaseValues[i] = m_ShData->GetSampleValues(databaseValues[i],m_SampleDirections);

    patientVectorValue = m_ShData->GetSampleValues(patientVectorValue,m_SampleDirections);

    return patientVectorValue.GetSize();
}

} // end namespace anima
