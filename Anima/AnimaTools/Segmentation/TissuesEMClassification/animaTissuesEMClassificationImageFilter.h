#pragma once

#include <iostream>
#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <vnl/vnl_matrix.h>
#include <vector>

namespace anima
{

template <typename TInputImage>
class TissuesEMClassificationImageFilter :
        public anima::MaskedImageToImageFilter< TInputImage, itk::VectorImage <double, 3> >
{
public:
    /** Standard class typedefs. */
    typedef TissuesEMClassificationImageFilter Self;
    typedef itk::VectorImage <double,3> TOutputImage;
    typedef anima::MaskedImageToImageFilter <TInputImage, TOutputImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(TissuesEMClassificationImageFilter, anima::MaskedImageToImageFilter)

    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef typename TInputImage::PixelType InputPixelType;

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::MaskImagePointer MaskImagePointer;

    /** My typedefs */
    typedef vnl_matrix <double> ParametersMatrixType;

    typedef itk::VectorImage <double,3> LocalPriorImageType;
    typedef LocalPriorImageType::Pointer LocalPriorImagePointer;

    typedef itk::Image <double,3> RealImageType;
    typedef RealImageType::Pointer RealImagePointer;

    /** Set/Get the mask on which to compute STAPLE estimate. */
    itkSetObjectMacro(LocalPriorImage, LocalPriorImageType)
    itkGetMacro(LocalPriorImage, LocalPriorImagePointer)

    /** Set/Get the maximum number of iterations after which the algorithm
      *  will be considered to have converged. */
    itkSetMacro(MaximumIterations, unsigned int)
    itkGetMacro(MaximumIterations, unsigned int)

    /** Set/Get the threshold for which a change in the maximum of the difference between parameters
      * shall be so small as to trigger termination of the estimation procedure.
      */
    itkSetMacro(RelativeConvergenceThreshold, double)
    itkGetMacro(RelativeConvergenceThreshold, double)

    /** Compute the M-step of the algorithm
      *  (i.e. the parameters of each expert from the reference standard and the data).
      */
    void EstimateClassesParameters();
    void EstimateClassesParameters(unsigned int classStart, unsigned int classEnd);

    bool endConditionReached();

    MaskImagePointer &GetClassificationAsLabelMap();
    RealImagePointer &GetZScoreMap();

    itkSetMacro(Verbose, bool)
    itkGetMacro(NumberOfClasses, unsigned int)

protected:
    TissuesEMClassificationImageFilter()
        : Superclass()
    {
        m_NumberOfClasses = 3;
        m_MaximumIterations = 100;
        m_Verbose = true;
    }

    virtual ~TissuesEMClassificationImageFilter() {}

    /** Compute the E-step of the algorithm
      *  (i.e. the reference standard from the parameters and the data).
      */
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    void GenerateOutputInformation() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void GenerateData() ITK_OVERRIDE;

    struct EMStepThreadStruct
    {
        Pointer Filter;
    };

    // Does the splitting and calls EstimatePerformanceParameters on a sub sample of experts
    static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreadEstimateClassesParams(void *arg);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(TissuesEMClassificationImageFilter);

    unsigned int m_MaximumIterations;
    double m_RelativeConvergenceThreshold;

    std::vector <double> m_ClassesPriors;
    unsigned int m_NumberOfClasses, m_NumberOfInputs;
    bool m_Verbose;

    std::vector < std::vector <double> > m_ClassesMeans, m_OldClassesMeans;
    std::vector < vnl_matrix <double> > m_ClassesVariances, m_InverseClassesVariances, m_OldClassesVariances;
    std::vector <double> m_ClassesVariancesSqrtDeterminants;

    LocalPriorImagePointer m_LocalPriorImage;
    MaskImagePointer m_LabelMap;
    RealImagePointer m_ZScoreMap;
};

} // end namespace anima

#include "animaTissuesEMClassificationImageFilter.hxx"
