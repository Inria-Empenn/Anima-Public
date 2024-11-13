#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{

template <class PixelScalarType>
class DTINonCentralChiEstimationImageFilter :
    public anima::MaskedImageToImageFilter< itk::Image<PixelScalarType,3>, itk::VectorImage<PixelScalarType,3> >
{
public:
    /** Standard class typedefs. */
    typedef DTINonCentralChiEstimationImageFilter<PixelScalarType> Self;
    typedef itk::Image<PixelScalarType,3> InputImageType;
    typedef itk::VectorImage<PixelScalarType,3> OutputImageType;
    typedef anima::MaskedImageToImageFilter< InputImageType, OutputImageType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(DTINonCentralChiEstimationImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::MaskImageType MaskImageType;

    itkSetMacro(EulerMascheroniConstant, double)

    void SetBValuesList(std::vector <double> bValuesList ) {m_BValuesList = bValuesList;}

    itkSetMacro(MaximumNumberOfBFGSIterations, unsigned int)
    itkSetMacro(MaximumNumberOfIterations, unsigned int)
    itkSetMacro(NumberOfCoils, unsigned int)

    itkSetMacro(StopThreshold, double)
    itkSetMacro(BFGSStopThreshold, double)
    itkSetMacro(B0Threshold, double)

    itkSetMacro(PValueThreshold, double)
    itkSetMacro(B0Index, unsigned int)

    itkSetMacro(InitialDTIImage, OutputImagePointer)
    itkSetMacro(InitialEstimatedB0Image, InputImagePointer)

    itkSetMacro(OptimizeB0Value, bool)
    itkSetMacro(RemoveDegeneratedTensors, bool)

    itkGetMacro(EffectiveCoilsImage, InputImageType *)
    itkGetMacro(LocalVarianceImage, InputImageType *)

    void AddGradientDirection(unsigned int i, std::vector <double> &grad);

protected:
    DTINonCentralChiEstimationImageFilter()
        : Superclass()
    {
        m_BValuesList.clear();
        m_GradientDirections.clear();

        m_EulerMascheroniConstant = 0.577;
        m_B0Threshold = 0;
        m_GlobalSigma = 0;
        m_BFGSStopThreshold = 10;
        m_StopThreshold = 1.0e-2;
        m_MaximumNumberOfBFGSIterations = 100;
        m_MaximumNumberOfIterations = 100;
        m_AverageBValue = 1000;
        m_NumberOfCoils = 1;

        m_OptimizeB0Value = false;

        m_InitialDTIImage = NULL;
        m_InitialEstimatedB0Image = NULL;

        m_EffectiveCoilsImage = NULL;
        m_LocalVarianceImage = NULL;
    }

    virtual ~DTINonCentralChiEstimationImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    void InitializeFromData(bool keepMask);
    double ComputeLocalSigma(unsigned int nbCoils, double oldLocalSigma, double b0Value, std::vector <double> &dtiValue, std::vector <double> &dwi);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(DTINonCentralChiEstimationImageFilter);

    double m_EulerMascheroniConstant;

    std::vector <double> m_BValuesList;
    double m_AverageBValue;
    std::vector< std::vector <double> > m_GradientDirections;

    OutputImagePointer m_InitialDTIImage;
    InputImagePointer m_InitialEstimatedB0Image;

    static const unsigned int m_NumberOfComponents = 6;

    //! Controls the number of coils, setting it to 1 goes back to Rician noise assumption
    unsigned int m_NumberOfCoils;

    unsigned int m_MaximumNumberOfBFGSIterations, m_MaximumNumberOfIterations;
    double m_BFGSStopThreshold;
    double m_StopThreshold;

    vnl_matrix <double> m_DesignMatrix;
    double m_GlobalSigma;
    bool m_OptimizeB0Value;

    // Global variance estimation parameters
    double m_PValueThreshold;
    double m_B0Threshold;
    unsigned int m_B0Index;

    bool m_RemoveDegeneratedTensors;

    InputImagePointer m_EffectiveCoilsImage;
    InputImagePointer m_LocalVarianceImage;
};

} // end namespace anima

#include "animaDTINonCentralChiEstimationImageFilter.hxx"
