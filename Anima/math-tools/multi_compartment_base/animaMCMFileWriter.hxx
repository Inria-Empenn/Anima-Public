#pragma once
#include "animaMCMFileWriter.h"

#include <itkImageRegionIterator.h>
#include <itkFileTools.h>
#include <animaReadWriteFunctions.h>

namespace anima
{

template <class PixelType, unsigned int ImageDimension>
MCMFileWriter <PixelType, ImageDimension>
::MCMFileWriter()
{
    m_FileName = "";
}

template <class PixelType, unsigned int ImageDimension>
MCMFileWriter <PixelType, ImageDimension>
::~MCMFileWriter()
{
}

template <class PixelType, unsigned int ImageDimension>
void
MCMFileWriter <PixelType, ImageDimension>
::SetFileName(std::string fileName)
{
    // if (fileName.find('.') != std::string::npos)
    //     fileName.erase(fileName.find_first_of('.'));
    m_FileName = fileName;
}

template <class PixelType, unsigned int ImageDimension>
void
MCMFileWriter <PixelType, ImageDimension>
::Update()
{
    if (m_FileName == "")
        throw itk::ExceptionObject(__FILE__, __LINE__,"No filename specified for writing",ITK_LOCATION);

    if (!m_InputImage)
        throw itk::ExceptionObject(__FILE__, __LINE__,"No input file to write",ITK_LOCATION);

    if (!m_InputImage->GetDescriptionModel())
        throw itk::ExceptionObject(__FILE__, __LINE__,"No reference model provided for writing MCM file",ITK_LOCATION);

    std::replace(m_FileName.begin(),m_FileName.end(),'\\','/');
    std::string noPathName = m_FileName;
    std::size_t lastSlashPos = m_FileName.find_last_of("/");

    if (lastSlashPos != std::string::npos)
    {
        noPathName.clear();
        noPathName.append(m_FileName.begin() + lastSlashPos + 1,m_FileName.end());
    }

    itk::FileTools::CreateDirectory(m_FileName.c_str());
    std::string headerName = m_FileName + ".mcm";
    std::ofstream outputHeaderFile(headerName.c_str());
    outputHeaderFile << "<?xml version=\"1.0\"?>" << std::endl;
    outputHeaderFile << "<Model>" << std::endl;

    // Output weights image
    ModelPointer descriptionModel = m_InputImage->GetDescriptionModel();

    unsigned int numberOfCompartments = descriptionModel->GetNumberOfCompartments();
    BaseOutputImagePointer weightsImage = BaseOutputImageType::New();
    weightsImage->Initialize();
    weightsImage->SetRegions(m_InputImage->GetLargestPossibleRegion());
    weightsImage->SetSpacing (m_InputImage->GetSpacing());
    weightsImage->SetOrigin (m_InputImage->GetOrigin());
    weightsImage->SetDirection (m_InputImage->GetDirection());
    weightsImage->SetVectorLength(numberOfCompartments);
    weightsImage->Allocate();

    typedef itk::ImageRegionIterator <InputImageType> InputImageIteratorType;
    typedef itk::ImageRegionIterator <BaseOutputImageType> BaseOutputImageIteratorType;
    typedef typename InputImageType::PixelType VectorType;

    BaseOutputImageIteratorType weightsItr(weightsImage,m_InputImage->GetLargestPossibleRegion());
    InputImageIteratorType inputItr(m_InputImage,m_InputImage->GetLargestPossibleRegion());

    VectorType tmpWeights(numberOfCompartments);
    VectorType workVector;
    while (!weightsItr.IsAtEnd())
    {
        workVector = inputItr.Get();

        for (unsigned int i = 0;i < numberOfCompartments;++i)
            tmpWeights[i] = workVector[i];

        weightsItr.Set(tmpWeights);

        ++weightsItr;
        ++inputItr;
    }

    std::string weightsName = m_FileName + "/";
    weightsName += noPathName;
    weightsName += "_weights.nrrd";

    anima::writeImage <BaseOutputImageType> (weightsName,weightsImage);
    std::string xmlFileNameWeights = noPathName + "_weights.nrrd";

    outputHeaderFile << "<Weights>" << xmlFileNameWeights << "</Weights>" << std::endl;

    // Start index is after weights for image data
    unsigned int pos = numberOfCompartments;

    for (unsigned int i = 0;i < descriptionModel->GetNumberOfCompartments();++i)
    {
        outputHeaderFile << "<Compartment>" << std::endl;
        switch(descriptionModel->GetCompartment(i)->GetCompartmentType())
        {
            case Stick:
                outputHeaderFile << "<Type>Stick</Type>" << std::endl;
                break;

            case Zeppelin:
                outputHeaderFile << "<Type>Zeppelin</Type>" << std::endl;
                break;

            case Tensor:
                outputHeaderFile << "<Type>Tensor</Type>" << std::endl;
                break;
                
            case NODDI:
                outputHeaderFile << "<Type>NODDI</Type>" << std::endl;
                break;

            case DDI:
                outputHeaderFile << "<Type>DDI</Type>" << std::endl;
                break;

            case FreeWater:
                outputHeaderFile << "<Type>FreeWater</Type>" << std::endl;
                break;

            case StationaryWater:
                outputHeaderFile << "<Type>StationaryWater</Type>" << std::endl;
                break;

            case Stanisz:
                outputHeaderFile << "<Type>Stanisz</Type>" << std::endl;
                break;

            case IsotropicRestrictedWater:
            default:
                outputHeaderFile << "<Type>IRWater</Type>" << std::endl;
                break;
        }

        // Output compartment image
        unsigned int compartmentSize = descriptionModel->GetCompartment(i)->GetCompartmentSize();

        BaseOutputImagePointer compartmentImage = BaseOutputImageType::New();
        compartmentImage->Initialize();
        compartmentImage->SetRegions(m_InputImage->GetLargestPossibleRegion());
        compartmentImage->SetSpacing (m_InputImage->GetSpacing());
        compartmentImage->SetOrigin (m_InputImage->GetOrigin());
        compartmentImage->SetDirection (m_InputImage->GetDirection());
        compartmentImage->SetVectorLength(compartmentSize);
        compartmentImage->Allocate();

        inputItr.GoToBegin();
        BaseOutputImageIteratorType compartmentItr(compartmentImage,m_InputImage->GetLargestPossibleRegion());

        VectorType tmpCompartment(compartmentSize);
        while (!compartmentItr.IsAtEnd())
        {
            workVector = inputItr.Get();

            for (unsigned int j = 0;j < compartmentSize;++j)
                tmpCompartment[j] = workVector[pos + j];

            compartmentItr.Set(tmpCompartment);

            ++compartmentItr;
            ++inputItr;
        }

        pos += compartmentSize;

        std::string fullPathCompartmentName = m_FileName + "/";
        std::string compartmentName = noPathName + "_";

        char tmpStr[2048];
        sprintf(tmpStr,"%d",i);
        compartmentName += tmpStr;
        compartmentName += ".nrrd";

        fullPathCompartmentName += compartmentName;

        std::string xmlFileNameCompartment = compartmentName;

        anima::writeImage <BaseOutputImageType> (fullPathCompartmentName,compartmentImage);
        outputHeaderFile << "<FileName>" << xmlFileNameCompartment << "</FileName>" << std::endl;

        outputHeaderFile << "</Compartment>" << std::endl;
    }

    outputHeaderFile << "</Model>" << std::endl;
    outputHeaderFile.close();
}

} // end namespace anima
