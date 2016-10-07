#include <tclap/CmdLine.h>

#include <animaResampleImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkAddImageFilter.h>
#include <itkImageRegionIterator.h>
#include <animaReadWriteFunctions.h>

int main(int ac, const char** av)
{
    std::string descriptionMessage = "Resampler tool for stitching several image into one.\n"
                                     "Applies linear transforms associated to the input images, then\n"
                                     "generates a larger image containing the mosaic of the transformed input images.\n"
                                     "INRIA / IRISA - VisAGeS Team";

    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);

    TCLAP::MultiArg<std::string> inArg("i","input","Input image (should be used several times)",true,"input image",cmd);
    TCLAP::MultiArg<std::string> trArg("t","trsf","Transformation (one linear transform for each input)",false,"transformation",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output mosaic image",true,"","output mosaic image",cmd);
    TCLAP::ValueArg<std::string> outMaskArg("O","out-mask","Output mosaic mask",false,"","output mosaic mask",cmd);
    TCLAP::ValueArg<std::string> geomArg("g","geometry","Geometry image",true,"","geometry image",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",
                                         false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // Find out the type of the image in file
    typedef itk::Image <float, 3> ImageType;
    ImageType::Pointer geometryImage = anima::readImage <ImageType> (geomArg.getValue());

    typedef itk::MatrixOffsetTransformBase <double,3> MatrixTransformType;
    typedef MatrixTransformType::Pointer MatrixTransformPointer;
    typedef itk::ContinuousIndex <double,3> ContinuousIndexType;
    ContinuousIndexType lowerCornerBoundingBox, upperCornerBoundingBox;

    for (unsigned int i = 0;i < ImageType::ImageDimension;++i)
    {
        lowerCornerBoundingBox[i] = geometryImage->GetLargestPossibleRegion().GetIndex()[i];
        upperCornerBoundingBox[i] = lowerCornerBoundingBox[i] + geometryImage->GetLargestPossibleRegion().GetSize()[i];
    }

    std::vector <std::string> inputImages = inArg.getValue();
    std::vector <std::string> inputTrsfs = trArg.getValue();

    // Explore images and transforms to get the boundaries of the future image
    for (unsigned int k = 0;k < inputImages.size();++k)
    {
        ImageType::Pointer input = anima::readImage <ImageType> (inputImages[k]);

        ImageType::IndexType lowerCornerIndex = input->GetLargestPossibleRegion().GetIndex();
        ImageType::IndexType upperCornerIndex = lowerCornerIndex + input->GetLargestPossibleRegion().GetSize();

        MatrixTransformPointer trsf = MatrixTransformType::New();
        trsf->SetIdentity();

        if (inputTrsfs.size() != 0)
        {
            itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
            reader->SetFileName(inputTrsfs[k]);
            reader->Update();

            const itk::TransformFileReader::TransformListType *trsfList = reader->GetTransformList();
            itk::TransformFileReader::TransformListType::const_iterator tr_it = trsfList->begin();

            trsf = dynamic_cast <MatrixTransformType *> ((*tr_it).GetPointer());
            MatrixTransformPointer tmpInvert = MatrixTransformType::New();
            trsf->GetInverse(tmpInvert);
            trsf = tmpInvert;
        }

        unsigned int numCorners = 1 << ImageType::ImageDimension;
        for (unsigned int i = 0;i < numCorners;++i)
        {
            ImageType::IndexType corner;
            unsigned int upper = i; // each bit indicates upper/lower on dim

            for(unsigned int dim = 0; dim < ImageType::ImageDimension; ++dim)
            {
                if (upper & 1)
                    corner[dim] = upperCornerIndex[dim];
                else
                    corner[dim] = lowerCornerIndex[dim];

                upper >>= 1;
            }

            ImageType::PointType tmpPoint;
            input->TransformIndexToPhysicalPoint(corner,tmpPoint);
            tmpPoint = trsf->TransformPoint(tmpPoint);

            itk::ContinuousIndex <double,3> tmpRes;
            geometryImage->TransformPhysicalPointToContinuousIndex(tmpPoint,tmpRes);
            for (unsigned int j = 0;j < ImageType::ImageDimension;++j)
            {
                if (tmpRes[j] < lowerCornerBoundingBox[j])
                    lowerCornerBoundingBox[j] = tmpRes[j];

                if (tmpRes[j] > upperCornerBoundingBox[j])
                    upperCornerBoundingBox[j] = tmpRes[j];
            }
        }
    }

    // Now define output geometry
    ImageType::Pointer outputImage = ImageType::New();
    outputImage->Initialize();

    for (unsigned int i = 0;i < ImageType::ImageDimension;++i)
    {
        lowerCornerBoundingBox[i] = std::floor(lowerCornerBoundingBox[i]);
        upperCornerBoundingBox[i] = std::ceil(upperCornerBoundingBox[i]);
    }

    ImageType::PointType origin;
    geometryImage->TransformContinuousIndexToPhysicalPoint(lowerCornerBoundingBox,origin);
    outputImage->SetOrigin(origin);
    outputImage->SetSpacing(geometryImage->GetSpacing());
    outputImage->SetDirection(geometryImage->GetDirection());

    ImageType::RegionType outputRegion;
    for (unsigned int i = 0;i < ImageType::ImageDimension;++i)
    {
        int lowerIndex = lowerCornerBoundingBox[i];
        int upperIndex = upperCornerBoundingBox[i];

        if (lowerIndex >= 0)
            outputRegion.SetIndex(i,lowerIndex);
        else
            outputRegion.SetIndex(i,0);

        outputRegion.SetSize(i,upperIndex - lowerIndex);
    }

    outputImage->SetRegions(outputRegion);
    outputImage->Allocate();
    outputImage->FillBuffer(0.0);

    typedef itk::Image <unsigned short,3> MaskImageType;
    MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetOrigin(origin);
    maskImage->SetSpacing(geometryImage->GetSpacing());
    maskImage->SetDirection(geometryImage->GetDirection());
    maskImage->SetRegions(outputRegion);
    maskImage->Allocate();
    maskImage->FillBuffer(0);

    // Now finally resample images in new geometry and average
    typedef anima::ResampleImageFilter <ImageType, ImageType> ResampleFilterType;
    typedef anima::ResampleImageFilter <MaskImageType, ImageType> MaskResampleFilterType;
    typedef itk::AddImageFilter <ImageType,ImageType,ImageType> AddFilterType;
    typedef itk::ImageRegionIterator <ImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskImageIteratorType;

    for (unsigned int i = 0;i < inputImages.size();++i)
    {
        std::cout << "Mosaicing image " << i+1 << ": " << inputImages[i];
        if (inputTrsfs.size() > 0)
            std::cout << " with transform " << inputTrsfs[i] << std::endl;
        else
            std::cout << std::endl;

        ImageType::Pointer input = anima::readImage <ImageType> (inputImages[i]);

        MatrixTransformPointer trsf = MatrixTransformType::New();
        trsf->SetIdentity();

        if (inputTrsfs.size() != 0)
        {
            itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
            reader->SetFileName(inputTrsfs[i]);
            reader->Update();

            const itk::TransformFileReader::TransformListType *trsfList = reader->GetTransformList();
            itk::TransformFileReader::TransformListType::const_iterator tr_it = trsfList->begin();

            trsf = dynamic_cast <MatrixTransformType *> ((*tr_it).GetPointer());
        }

        ResampleFilterType::Pointer imageResampler = ResampleFilterType::New();
        imageResampler->SetTransform(trsf);
        imageResampler->SetSize(outputRegion.GetSize());
        imageResampler->SetOutputOrigin(origin);
        imageResampler->SetOutputSpacing(outputImage->GetSpacing());
        imageResampler->SetOutputDirection(outputImage->GetDirection());
        imageResampler->SetInput(input);
        imageResampler->SetNumberOfThreads(nbpArg.getValue());
        imageResampler->Update();

        AddFilterType::Pointer addFilter = AddFilterType::New();
        addFilter->SetInput1(outputImage);
        addFilter->SetInput2(imageResampler->GetOutput());
        addFilter->SetNumberOfThreads(nbpArg.getValue());
        addFilter->Update();

        outputImage = addFilter->GetOutput();
        outputImage->DisconnectPipeline();

        MaskImageType::Pointer tmpMask = MaskImageType::New();
        tmpMask->Initialize();
        tmpMask->SetOrigin(input->GetOrigin());
        tmpMask->SetSpacing(input->GetSpacing());
        tmpMask->SetDirection(input->GetDirection());
        tmpMask->SetRegions(input->GetLargestPossibleRegion());
        tmpMask->Allocate();
        tmpMask->FillBuffer(1);

        MaskResampleFilterType::Pointer maskImResampler = MaskResampleFilterType::New();
        maskImResampler->SetTransform(trsf);
        maskImResampler->SetSize(outputRegion.GetSize());
        maskImResampler->SetOutputOrigin(origin);
        maskImResampler->SetOutputSpacing(outputImage->GetSpacing());
        maskImResampler->SetOutputDirection(outputImage->GetDirection());
        maskImResampler->SetInput(tmpMask);
        maskImResampler->SetNumberOfThreads(nbpArg.getValue());
        maskImResampler->Update();

        MaskImageIteratorType outputMaskItr(maskImage,maskImage->GetLargestPossibleRegion());
        ImageIteratorType trsfMaskItr(maskImResampler->GetOutput(),maskImage->GetLargestPossibleRegion());

        while (!outputMaskItr.IsAtEnd())
        {
            if (trsfMaskItr.Get() > 0)
                outputMaskItr.Set(outputMaskItr.Get() + 1);

            ++trsfMaskItr;
            ++outputMaskItr;
        }
    }

    // Finally divide by the number of images at each pixel
    ImageIteratorType outItr(outputImage,outputImage->GetLargestPossibleRegion());
    MaskImageIteratorType maskItr(maskImage,maskImage->GetLargestPossibleRegion());
    while (!outItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
            outItr.Set(outItr.Get() / maskItr.Get());
        else
            outItr.Set(0);

        ++outItr;
        ++maskItr;
    }

    anima::writeImage <ImageType> (outArg.getValue(),outputImage);

    if (outMaskArg.getValue() != "")
        anima::writeImage <MaskImageType> (outMaskArg.getValue(),maskImage);

    return EXIT_SUCCESS;
}
