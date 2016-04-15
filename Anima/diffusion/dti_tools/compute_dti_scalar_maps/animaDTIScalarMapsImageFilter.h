#pragma once

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include <itkConfigure.h>

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima

{
/** \class DTIScalarMapsImageFilter
 * \brief Applies an variance filter to an image
 *
 * Computes tow image where a given pixel is the FA or ADC value of the
 * tensor of the corresponding input pixel.
 *
 */
template <unsigned int ImageDimension = 3>
class DTIScalarMapsImageFilter :
    public itk::ImageToImageFilter< itk::VectorImage <float, ImageDimension>,
                                    itk::Image <float, ImageDimension> >
{
public:

    /** Convenient typedefs for simplifying declarations. */
    typedef itk::VectorImage <float, ImageDimension> InputImageType;
    typedef itk::Image< float, ImageDimension> OutputImageType;
    typedef InputImageType TensorImageType;
    typedef OutputImageType FAImageType;
    typedef OutputImageType ADCImageType;

    /** Extract dimension from input and output image. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        InputImageType::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        OutputImageType::ImageDimension);


    /** Standard class typedefs. */
    typedef DTIScalarMapsImageFilter                                      Self;
    typedef itk::ImageToImageFilter< InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self>                                   Pointer;
    typedef itk::SmartPointer<const Self>                             ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(DTIScalarMapsImageFilter, ImageToImageFilter)

    /** Image typedef support. */
    typedef typename TensorImageType::PixelType               TensorVectorType;

    typedef typename TensorImageType::RegionType  TensorImageRegionType;
    typedef typename OutputImageType::RegionType  OutputImageRegionType;
    typedef typename TensorImageType::SizeType    TensorImageSizeType;

    /**  Create the Output */
    itk::DataObject::Pointer MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

    typename OutputImageType::Pointer GetADCImage() {return this->GetOutput(0);}
    typename OutputImageType::Pointer GetFAImage() {return this->GetOutput(1);}
    typename OutputImageType::Pointer GetAxialDiffusivityImage() {return this->GetOutput(2);}
    typename OutputImageType::Pointer GetRadialDiffusivityImage() {return this->GetOutput(3);}

protected:
    DTIScalarMapsImageFilter();
    virtual ~DTIScalarMapsImageFilter() {}

    /** DTIScalarMapsImageFilter can be implemented as a multithreaded filter.
     * Therefore, this implementation provides a ThreadedGenerateData()
     * routine which is called for each processing thread. The output
     * image data is allocated automatically by the superclass prior to
     * calling ThreadedGenerateData().  ThreadedGenerateData can only
     * write to the portion of the output image specified by the
     * parameter "outputRegionForThread"
     *
     * \sa ImageToImageFilter::ThreadedGenerateData(),
     *     ImageToImageFilter::GenerateData() */
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                              itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    DTIScalarMapsImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented


};

} // end of namespace anima

#include "animaDTIScalarMapsImageFilter.hxx"
