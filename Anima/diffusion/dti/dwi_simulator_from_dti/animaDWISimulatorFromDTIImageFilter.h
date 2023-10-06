#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{

template <class PixelScalarType>
class DWISimulatorFromDTIImageFilter :
        public itk::ImageToImageFilter < itk::VectorImage<PixelScalarType,3>, itk::VectorImage<PixelScalarType,3> >
{
public:
    /** Standard class typedefs. */
    typedef DWISimulatorFromDTIImageFilter<PixelScalarType> Self;
    typedef itk::VectorImage<PixelScalarType,3> TInputImage;
    typedef TInputImage OutputB0ImageType;
    typedef itk::VectorImage<PixelScalarType,3> DTIImageType;
    typedef itk::Image<PixelScalarType,4> Image4DType;
    typedef itk::Image<PixelScalarType,3> S0ImageType;
    typedef typename S0ImageType::Pointer S0ImagePointer;
    typedef itk::VectorImage<PixelScalarType,3> TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(DWISimulatorFromDTIImageFilter, ImageToImageFilter)

    /** Image typedef support */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::PixelType InputPixelType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetBValuesList (std::vector <double> bValuesList) {m_BValuesList = bValuesList;}

    void AddGradientDirection(unsigned int i, std::vector <double> &grad);

    itkGetMacro(S0Value, double)
    itkSetMacro(S0Value, double)

    itkSetMacro(S0Image, S0ImagePointer)

    Image4DType *GetOutputAs4DImage();

protected:
    DWISimulatorFromDTIImageFilter()
        : Superclass()
    {
        m_BValuesList.clear();
        m_GradientDirections.clear();

        m_S0Value = 200;
        m_S0Image = 0;
    }

    virtual ~DWISimulatorFromDTIImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    bool isZero(InputPixelType &vec);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(DWISimulatorFromDTIImageFilter);

    std::vector <double> m_BValuesList;
    std::vector < std::vector <double> > m_GradientDirections;

    double m_S0Value;
    S0ImagePointer m_S0Image;

    static const unsigned int m_NumberOfComponents = 6;

    typename Image4DType::Pointer m_Output4D;
};

} // end namespace anima

#include "animaDWISimulatorFromDTIImageFilter.hxx"
