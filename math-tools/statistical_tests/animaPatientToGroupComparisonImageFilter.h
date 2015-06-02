#pragma once

#include <iostream>
#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <vector>

namespace anima
{

/**
     * \brief Implements patient to group comparison as in http://dx.doi.org/10.1007/978-3-540-85988-8_116
     *
     * Provides an ITK implementation of z-score and p-value computation telling how much different a single vector image is
     * from a database of controls. It can virtually process any type of input, however for tensor images, be warned
     * that they should be expressed as log-vectors.
     *
     */
template <class PixelScalarType>
class PatientToGroupComparisonImageFilter :
        public anima::MaskedImageToImageFilter< itk::VectorImage <PixelScalarType, 3> , itk::Image <PixelScalarType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef PatientToGroupComparisonImageFilter<PixelScalarType> Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(PatientToGroupComparisonImageFilter, MaskedImageToImageFilter);

    /** Image typedef support */
    typedef itk::VectorImage <PixelScalarType, 3> InputImageType;
    typedef itk::Image <PixelScalarType, 3> OutputImageType;

    typedef typename InputImageType::PixelType VectorType;

    enum TestType
    {
        CHI_SQUARE = 0,
        FISHER
    };

    typedef anima::MaskedImageToImageFilter< InputImageType, OutputImageType > Superclass;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    typedef vnl_matrix<double> CovMatrixType;

    void AddDatabaseInput(InputImageType *tmpIm)
    {
        m_DatabaseImages.push_back(tmpIm);
    }

    itkSetMacro(NumEigenValuesPCA, unsigned int);
    itkSetMacro(ExplainedRatio, double);

    itkSetMacro(StatisticalTestType, TestType);

protected:
    PatientToGroupComparisonImageFilter()
        : Superclass()
    {
        this->SetNumberOfRequiredOutputs(2);
        this->SetNthOutput(0,this->MakeOutput(0));
        this->SetNthOutput(1,this->MakeOutput(1));

        m_ExplainedRatio = 0.9;
        m_NumEigenValuesPCA = 6;

        m_DatabaseImages.clear();
        m_StatisticalTestType = FISHER;
    }

    virtual ~PatientToGroupComparisonImageFilter() {}

    virtual void BeforeThreadedGenerateData(void);
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

    virtual unsigned int SampleFromDiffusionModels(std::vector <VectorType> &databaseValues, VectorType &patientVectorValue)
    {
        // Reimplement in subclasses if they need to sample from the distribution before doing the projection and comparison
        return patientVectorValue.GetSize();
    }

private:
    PatientToGroupComparisonImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    bool isZero(const VectorType &vec);

    unsigned int GetPCAVectorsFromData(std::vector < itk::VariableLengthVector <double> > &databaseVectors,
                                       itk::VariableLengthVector <double> &patientVectorValue);

    unsigned int m_NumEigenValuesPCA;
    double m_ExplainedRatio;

    std::vector <InputImagePointer> m_DatabaseImages;
    TestType m_StatisticalTestType;
};

} // end namespace anima

#include "animaPatientToGroupComparisonImageFilter.hxx"
