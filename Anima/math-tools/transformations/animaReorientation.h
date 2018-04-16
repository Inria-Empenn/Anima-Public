#pragma once

#include <itkImage.h>
#include <itkOrientImageFilter.h>
#include <animaGradientFileReader.h>
#include <itkExtractImageFilter.h>

namespace anima
{
template <class ImageType>
typename itk::SmartPointer<ImageType>
reorient3DImage(typename itk::SmartPointer<ImageType> input,
              itk::SpatialOrientation::ValidCoordinateOrientationFlags orientation)
{
    // Type redefinition and cast are only here for compilation purpose.
    typedef itk::Image<typename ImageType::PixelType, 3> CheckedImageType;
    typename CheckedImageType::Pointer checkedInput;
    checkedInput = dynamic_cast<CheckedImageType *>(input.GetPointer());

    typedef itk::OrientImageFilter<CheckedImageType, CheckedImageType> OrientImageFilter;
    typename OrientImageFilter::Pointer orienter = OrientImageFilter::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(orientation);
    orienter->SetInput(checkedInput);
    orienter->Update();
    return dynamic_cast<ImageType *>(orienter->GetOutput());
}

template <class ImageType>
typename itk::SmartPointer<ImageType>
reorientImage(typename itk::SmartPointer<ImageType> input,
              itk::SpatialOrientation::ValidCoordinateOrientationFlags orientation)
{
    if(ImageType::ImageDimension == 4)
    {
        typedef itk::Image<typename ImageType::PixelType, 3> ImageToReorientType;
        typedef itk::ExtractImageFilter<ImageType, ImageToReorientType> ExtractImageFilterType;
        typename ExtractImageFilterType::Pointer extractFilter = ExtractImageFilterType::New();
        extractFilter->SetInput(input);
        extractFilter->SetDirectionCollapseToGuess();

        typename ImageType::IndexType extractIndex;
        typename ImageType::SizeType extractSize;
        typename ImageType::RegionType extractRegion;

        for (unsigned int d = 0; d < 3; ++d)
        {
            extractIndex[d] = 0;
            extractSize[d] = input->GetLargestPossibleRegion().GetSize()[d];
        }
        extractIndex[3] = 0;
        extractRegion.SetIndex(extractIndex);
        extractSize[3] = 0;
        extractRegion.SetSize(extractSize);
        extractFilter->SetExtractionRegion(extractRegion);

        extractFilter->Update();

        typename ImageToReorientType::Pointer reorientedBase;
        reorientedBase = reorient3DImage<ImageToReorientType>(extractFilter->GetOutput(), orientation);
        typename ImageType::PointType reorientedOrigin;
        typename ImageType::SpacingType reorientedSpacing;
        typename ImageType::DirectionType reorientedDirection;
        typename ImageType::IndexType reorientedStart;
        typename ImageType::SizeType reorientedEnd;


        for(unsigned int d = 0; d < 3; ++d)
        {
            reorientedOrigin[d] = reorientedBase->GetOrigin()[d];
            reorientedSpacing[d] = reorientedBase->GetSpacing()[d];
            for(unsigned int e = 0; e < 3; ++e)
                reorientedDirection[d][e] = reorientedBase->GetDirection()[d][e];
            reorientedDirection[d][3] = input->GetDirection()[d][3];

            reorientedStart[d] = reorientedBase->GetLargestPossibleRegion().GetIndex()[d];
            reorientedEnd[d] = reorientedBase->GetLargestPossibleRegion().GetSize()[d];
        }
        reorientedOrigin[3] = input->GetOrigin()[3];
        reorientedSpacing[3] = input->GetSpacing()[3];
        for(unsigned int f = 0; f < 4; ++f)
            reorientedDirection[3][f] = input->GetDirection()[3][f];
        reorientedStart[3] = input->GetLargestPossibleRegion().GetIndex()[3];
        reorientedEnd[3] = input->GetLargestPossibleRegion().GetSize()[3];

        typename ImageType::Pointer reorientedImage = ImageType::New();
        reorientedImage->SetRegions(typename ImageType::RegionType(reorientedStart, reorientedEnd));
        reorientedImage->SetOrigin(reorientedOrigin);
        reorientedImage->SetSpacing(reorientedSpacing);
        reorientedImage->SetDirection(reorientedDirection);
        reorientedImage->Allocate();

        typedef itk::ImageRegionIterator <ImageType> FillIteratorType;
        FillIteratorType fillItr(reorientedImage, reorientedImage->GetLargestPossibleRegion());

        typedef itk::ImageRegionConstIterator <ImageToReorientType> ReorientedIteratorType;
        ReorientedIteratorType reoItr(reorientedBase, reorientedBase->GetLargestPossibleRegion());
        while (!reoItr.IsAtEnd())
        {
            fillItr.Set(reoItr.Get());
            ++fillItr;
            ++reoItr;
        }

        for(unsigned int i = 1; i < input->GetLargestPossibleRegion().GetSize()[3]; ++i)
        {
            extractIndex[3] = i;
            extractRegion.SetIndex(extractIndex);
            extractFilter->SetExtractionRegion(extractRegion);
            extractFilter->SetDirectionCollapseToGuess();
            extractFilter->Update();

            reorientedBase = reorient3DImage<ImageToReorientType>(extractFilter->GetOutput(), orientation);
            reoItr = ReorientedIteratorType(reorientedBase, reorientedBase->GetLargestPossibleRegion());
            while (!reoItr.IsAtEnd())
            {
                fillItr.Set(reoItr.Get());
                ++fillItr;
                ++reoItr;
            }
        }
        return reorientedImage;
    }
    else if(ImageType::ImageDimension == 3)
    {
        return reorient3DImage<ImageType>(input, orientation);
    }
    else
    {
        unsigned int dimension = ImageType::ImageDimension;
        std::string msg  = "Reorient an image of dimension "
                + std::to_string(dimension)
                + "is not supported yet.";
        throw itk::ExceptionObject(__FILE__, __LINE__, msg , ITK_LOCATION);
    }
}

template <class ImageType, class GradientType>
void reorientGradients(typename itk::SmartPointer<ImageType> input, std::vector<GradientType> &gradients)
{
    // What we want to do here is reorient a list of gradient
    // vectors in order to have them in image coordinates from a set of gradients in real coordinates

    // To do so we apply 'reorientedGrad = D^-1 o grad'

    typename ImageType::DirectionType inversedD1 = input->GetInverseDirection();

    for(unsigned int g = 0; g < gradients.size(); g++)
    {
        GradientType reorientedGrad(0.0);
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
                reorientedGrad[i] += inversedD1[i][j] * gradients[g][j];
        }
        gradients[g] = reorientedGrad;
    }
}

} // end namespace anima

