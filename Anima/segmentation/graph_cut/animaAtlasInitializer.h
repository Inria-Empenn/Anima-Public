#pragma once

#include "animaModelInitializer.h"
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

/** Class initializing 3 gaussian distributions from an atlas
   *
   */
template <typename TInputImage, typename TMaskImage, typename TAtlasImage = TInputImage>
class AtlasInitializer: public ModelInitializer
{
public:

    /** Standard class typedefs. */
    typedef AtlasInitializer  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(AtlasInitializer, itk::ProcessObject);

    /**  Type of the input image. */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< InputImageType > InputConstIteratorType;

    /**  Type of the mask image. */
    typedef TMaskImage MaskImageType;
    typedef typename MaskImageType::ConstPointer MaskImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< MaskImageType > MaskConstIteratorType;


    /**  Type of the mask image. */
    typedef TAtlasImage AtlasImageType;
    typedef typename AtlasImageType::ConstPointer AtlasImageConstPointer;
    typedef typename itk::ImageRegionIterator< AtlasImageType > AtlasIteratorType;
    typedef typename itk::ImageRegionConstIterator< AtlasImageType > AtlasConstIteratorType;

    typedef double                    NumericType;
    typedef itk::VariableSizeMatrix<NumericType> DoubleVariableSizeMatrixType;


    /** @brief calculate parameters of gaussians from images
     *  with the probability of the different classes [csf] [gm] [wm]
      */
    void Update() ITK_OVERRIDE;

    /** The mri images.*/
    void SetInputImage1(const TInputImage* image);
    void SetInputImage2(const TInputImage* image);
    void SetInputImage3(const TInputImage* image);

    /** The atlas images.*/
    void SetAtlasImage1(const TAtlasImage* image);
    void SetAtlasImage2(const TAtlasImage* image);
    void SetAtlasImage3(const TAtlasImage* image);

    /** mask in which the segmentation will be performed
       */
    void SetMask(const TMaskImage* mask);
    typename TMaskImage::ConstPointer GetMask();


protected:

    AtlasInitializer(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    typename TAtlasImage::ConstPointer GetAtlasImage1();
    typename TAtlasImage::ConstPointer GetAtlasImage2();
    typename TAtlasImage::ConstPointer GetAtlasImage3();

    typename TInputImage::ConstPointer GetInputImage1();
    typename TInputImage::ConstPointer GetInputImage2();
    typename TInputImage::ConstPointer GetInputImage3();

    AtlasInitializer()
    {
        this->SetNumberOfRequiredInputs(7);
        m_Modalities = 3;
    }
    virtual ~AtlasInitializer(){}

    std::vector<InputImageConstPointer > m_ImagesVector;
    std::vector<AtlasImageConstPointer > m_AtlasVector;
    unsigned int m_Modalities;

};

}

#include "animaAtlasInitializer.hxx"
