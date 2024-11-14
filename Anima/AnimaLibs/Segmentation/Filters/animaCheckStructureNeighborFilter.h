#pragma once

#include <animaReadWriteFunctions.h>

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelContourImageFilter.h>

namespace anima
{
/** @brief Class removing lesions that are not sufficiently in the white matter.
 * Intensity rules may not be enough to discard false positives, therefore we also use localization information.
 * Considering that MS lesions are typically located in WM, we remove the detected ones that do not sufficiently achieve this condition.
 * This filter take two entries: a lesion mask and a white matter map.
 *
 */
template<typename TInput, typename TMask, typename TOutput = TInput>
class CheckStructureNeighborFilter :
public itk::ImageToImageFilter< TInput , TOutput >
{
public:
    /** Standard class typedefs. */
    typedef CheckStructureNeighborFilter Self;
    typedef itk::ImageToImageFilter< TInput , TOutput > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(CheckStructureNeighborFilter, ImageToImageFilter)

    /** Image typedef support */

    /**  Type of the input image. */
    typedef TInput InputImageType;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename itk::ImageRegionConstIterator< InputImageType > InputConstIteratorType;

    /**  Type of the map image. */
    typedef TMask MaskImageType;
    typedef typename MaskImageType::ConstPointer MaskImageConstPointer;
    typedef typename MaskImageType::PixelType MaskPixelType;
    typedef typename itk::ImageRegionConstIterator< MaskImageType > MaskConstIteratorType;

    /**  Type of the mask image. */
    typedef TOutput OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename OutputImageType::PixelType OutputPixelType;
    typedef typename itk::ImageRegionIterator <OutputImageType> OutputIteratorType;

    typedef int PixelTypeInt;
    typedef itk::Image <PixelTypeInt,3> ImageTypeInt;
    typedef itk::ImageRegionIterator <ImageTypeInt> ImageIteratorTypeInt;

    typedef typename itk::ConnectedComponentImageFilter <InputImageType,ImageTypeInt> ConnectedComponentFilterType;
    typedef itk::LabelContourImageFilter<ImageTypeInt,ImageTypeInt> LabelContourFilterType;
    typedef itk::BinaryBallStructuringElement<PixelTypeInt, 3> StructuringElementType;
    typedef itk::GrayscaleDilateImageFilter <ImageTypeInt,ImageTypeInt,StructuringElementType> DilateFilterType;

    /** The mri images.*/
    void SetInputClassification(const TInput* image);
    void SetInputMap(const TMask* image);

    void WriteOutputs();

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
        m_Tol = tol;
    }

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(LabelToCheck, int)
    itkGetMacro(LabelToCheck, int)

    itkSetMacro(Verbose, bool)
    itkGetMacro(Verbose, bool)

    itkSetMacro(Ratio, double)
    itkGetMacro(Ratio, double)

    std::string GetOutputFilename() {return m_OutputFilename;}
    void SetOutputFilename(std::string filename) {m_OutputFilename=filename;}

protected:
    CheckStructureNeighborFilter()
    {
        this->SetNumberOfRequiredOutputs(1);
        this->SetNumberOfRequiredInputs(2);

        m_Verbose=false;
        m_Tol = 0.0001;
        m_Ratio = 0;

        this->SetNumberOfWorkUnits(itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());
    }

    virtual ~CheckStructureNeighborFilter()
    {
    }

    typename TInput::ConstPointer GetInputClassification();
    typename TMask::ConstPointer GetInputMap();

    void GenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(CheckStructureNeighborFilter);

    std::string m_OutputFilename;
    int m_LabelToCheck;
    bool m_Verbose;
    double m_Tol;
    double m_Ratio;

};

} // end of namespace anima

#include "animaCheckStructureNeighborFilter.hxx"
