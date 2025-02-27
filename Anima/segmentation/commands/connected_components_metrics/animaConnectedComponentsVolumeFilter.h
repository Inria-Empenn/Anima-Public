#pragma once

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <animaReadWriteFunctions.h>

namespace anima
{

template <typename ImageType>
class ConnectedComponentsVolumeFilter :
public itk::ImageToImageFilter< ImageType , ImageType >
{
public:
    /** Standard class typedefs. */
    typedef ConnectedComponentsVolumeFilter Self;
    typedef itk::ImageToImageFilter< ImageType, ImageType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(ConnectedComponentsVolumeFilter, ImageToImageFilter)

    /** Image typedef support */

    /**  Type of image. */
    typedef typename ImageType::Pointer ImagePointerType;
    typedef typename ImageType::ConstPointer ImageConstPointerType;
    typedef typename ImageType::PixelType ImagePixelType;
    typedef typename ImageType::IndexType ImageIndexType;
    typedef typename itk::ImageRegionIterator< ImageType > ImageIteratorType;
    typedef typename itk::ImageRegionConstIterator< ImageType > ImageConstIteratorType;

    typedef typename itk::ConnectedComponentImageFilter <ImageType,ImageType> ConnectedComponentFilterType;
    typedef typename itk::RelabelComponentImageFilter <ImageType,ImageType> RelabelComponentFilterType;

    typedef typename ImageType::SpacingType spacingType;
    typedef typename ImageType::SpacingValueType spacingValueType;

    void WriteOutputs();

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(MinSizeMM3, double)
    itkGetMacro(MinSizeMM3, double)

    itkSetMacro(FullyConnected, bool)
    itkGetMacro(FullyConnected, bool)

    itkSetMacro(Verbose, bool)
    itkGetMacro(Verbose, bool)

    itkSetMacro(Dimension, unsigned int)
    itkGetMacro(Dimension, unsigned int)

    itkGetMacro(OriginalNumberOfObjects, unsigned int)
    itkGetMacro(NumberOfObjects, unsigned int)
    itkGetMacro(TotalVolume, double)

    void SetOutputFilename(std::string outputFilename) {m_Outputfilename=outputFilename;}
    std::string GetOutputFilename() {return m_Outputfilename;}

    void SetOutputVolumeFilename(std::string outputFilename) {m_OutputVolumefilename=outputFilename;}
    std::string GetOutputVolumefFilename() {return m_OutputVolumefilename;}

    std::vector<double> GetVolumeVector() {return m_vect_volume_lesionMM3;}

    spacingType GetSpacing(){return m_Spacing;}
    spacingValueType GetSpacingTot(){return m_SpacingTot;}

protected:
    ConnectedComponentsVolumeFilter()
    {
        this->SetNumberOfRequiredOutputs(1);
        this->SetNumberOfRequiredInputs(1);

        m_MinSizeMM3 = 0;
        m_MinSizeVoxel = 0;
        m_Dimension = 2;
        m_FullyConnected = false;
        m_Verbose = false;

        m_OriginalNumberOfObjects = 0;
        m_NumberOfObjects = 0;
        m_TotalVolume= 0;
    }

    virtual ~ConnectedComponentsVolumeFilter()
    {
    }

    void Print(std::ostream& os);
    void ComputeSpacing();
    void ComputeMinimumLesionSize();
    void GenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(ConnectedComponentsVolumeFilter);

    // Parameters
    double m_MinSizeMM3;
    unsigned int m_MinSizeVoxel;
    unsigned int m_Dimension;
    bool m_FullyConnected;
    bool m_Verbose;

    // Output
    std::string m_Outputfilename;
    std::string m_OutputVolumefilename;
    unsigned int m_OriginalNumberOfObjects;
    unsigned int m_NumberOfObjects;
    double m_TotalVolume;
    std::vector<double> m_vect_volume_lesionMM3;
    spacingType m_Spacing;
    spacingValueType m_SpacingTot;
};

} // end of namespace anima

#include "animaConnectedComponentsVolumeFilter.hxx"
