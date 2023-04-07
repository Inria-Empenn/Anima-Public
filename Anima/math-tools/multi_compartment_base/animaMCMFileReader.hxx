#pragma once
#include "animaMCMFileReader.h"

#include <animaReadWriteFunctions.h>

#include <animaFreeWaterCompartment.h>
#include <animaIsotropicRestrictedWaterCompartment.h>
#include <animaStationaryWaterCompartment.h>
#include <animaStaniszCompartment.h>
#include <animaStickCompartment.h>
#include <animaZeppelinCompartment.h>
#include <animaTensorCompartment.h>
#include <animaNODDICompartment.h>

#include <itkImageRegionIterator.h>
#include <tinyxml2.h>

namespace anima
{

template <class PixelType, unsigned int ImageDimension>
MCMFileReader <PixelType, ImageDimension>
::MCMFileReader()
{
    m_FileName = "";
}

template <class PixelType, unsigned int ImageDimension>
MCMFileReader <PixelType, ImageDimension>
::~MCMFileReader()
{
}

template <class PixelType, unsigned int ImageDimension>
void
MCMFileReader <PixelType, ImageDimension>
::Update()
{
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLError loadOk = doc.LoadFile(m_FileName.c_str());

    std::replace(m_FileName.begin(),m_FileName.end(),'\\','/');
    std::string basePath, baseName;
    std::size_t lastSlashPos = m_FileName.find_last_of("/");

    if (lastSlashPos != std::string::npos)
    {
        basePath.append(m_FileName.begin(),m_FileName.begin() + lastSlashPos);
        ++lastSlashPos;
    }
    else
        lastSlashPos = 0;

    if (basePath.length() != 0)
        basePath = basePath + "/";

    std::size_t lastDotPos = m_FileName.find_last_of(".");
    baseName.append(m_FileName.begin() + lastSlashPos, m_FileName.begin() + lastDotPos);

    if (loadOk != tinyxml2::XML_SUCCESS)
    {
        std::string error("Unable to read input summary file: ");
        error += m_FileName;
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    // Read XML file into transform information vector
    tinyxml2::XMLElement *modelNode = doc.FirstChildElement( "Model" );
    if (!modelNode)
    {
        std::string error("Model node missing in ");
        error += m_FileName;
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    tinyxml2::XMLElement *weightsNode = modelNode->FirstChildElement( "Weights" );

    if (!weightsNode)
    {
        std::string error("Weights file missing in ");
        error += m_FileName;
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }
    
    std::string weightsFileName = basePath + baseName + "/";
    weightsFileName += weightsNode->GetText();
    
    BaseInputImagePointer weightsImage = anima::readImage <BaseInputImageType>(weightsFileName);
    unsigned int numCompartments = weightsImage->GetVectorLength();
    
    tinyxml2::XMLElement *compartmentNode = modelNode->FirstChildElement( "Compartment" );
    ModelPointer referenceModel = ModelType::New();
    
    std::vector <BaseInputImagePointer> compartmentImages;
    while (compartmentNode)
    {
        tinyxml2::XMLElement *typeNode = compartmentNode->FirstChildElement("Type");
        std::string compartmentType = typeNode->GetText();
        
        anima::BaseCompartment::Pointer additionalCompartment = this->CreateCompartmentForType(compartmentType);
        
        referenceModel->AddCompartment(1.0 / numCompartments,additionalCompartment);
        
        tinyxml2::XMLElement *fileNameNode = compartmentNode->FirstChildElement("FileName");
        std::string imageFileName = basePath + baseName + "/";
        imageFileName += fileNameNode->GetText();

        compartmentImages.push_back(anima::readImage <BaseInputImageType> (imageFileName));

        compartmentNode = compartmentNode->NextSiblingElement("Compartment");
    }
    
    unsigned int vectorFinalSize = referenceModel->GetSize();
    m_OutputImage = OutputImageType::New();
    m_OutputImage->Initialize();
    m_OutputImage->CopyInformation(compartmentImages[0]);
    m_OutputImage->SetRegions(compartmentImages[0]->GetLargestPossibleRegion());
    m_OutputImage->SetNumberOfComponentsPerPixel(vectorFinalSize);
    m_OutputImage->Allocate();
    m_OutputImage->SetDescriptionModel(referenceModel);
    
    typedef itk::ImageRegionIterator <BaseInputImageType> InputImageIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutputImageIteratorType;

    std::vector <InputImageIteratorType> inputIterators(numCompartments);
    for (unsigned int i = 0;i < numCompartments;++i)
        inputIterators[i] = InputImageIteratorType(compartmentImages[i],m_OutputImage->GetLargestPossibleRegion());
    
    InputImageIteratorType weightsIterator(weightsImage,m_OutputImage->GetLargestPossibleRegion());
    OutputImageIteratorType outItr(m_OutputImage,m_OutputImage->GetLargestPossibleRegion());

    typedef typename OutputImageType::PixelType ImagePixelType;
    ImagePixelType outValue(vectorFinalSize);
    ImagePixelType weightsData(vectorFinalSize);

    ImagePixelType inputCompartmentValue;
    ModelType::ModelOutputVectorType inputMCMValue;
    ModelType::ListType weightsVector(numCompartments);

    while (!outItr.IsAtEnd())
    {
        weightsData = weightsIterator.Get();
        for (unsigned int i = 0;i < numCompartments;++i)
            weightsVector[i] = weightsData[i];

        referenceModel->SetCompartmentWeights(weightsVector);

        for (unsigned int i = 0;i < numCompartments;++i)
        {
            inputCompartmentValue = inputIterators[i].Get();
            unsigned int compartmentSize = inputCompartmentValue.GetSize();

            if (inputMCMValue.GetSize() != compartmentSize)
                inputMCMValue.SetSize(compartmentSize);
            for (unsigned int j = 0;j < compartmentSize;++j)
                inputMCMValue[j] = inputCompartmentValue[j];

            referenceModel->GetCompartment(i)->SetCompartmentVector(inputMCMValue);
        }

        outValue = referenceModel->GetModelVector();
        outItr.Set(outValue);

        ++outItr;
        ++weightsIterator;
        for (unsigned int i = 0;i < numCompartments;++i)
            ++inputIterators[i];
    }
}

template <class PixelType, unsigned int ImageDimension>
anima::BaseCompartment::Pointer
MCMFileReader <PixelType, ImageDimension>
::CreateCompartmentForType(std::string &compartmentType)
{
    anima::BaseCompartment::Pointer additionalCompartment;

    if (compartmentType == "FreeWater")
        additionalCompartment = anima::FreeWaterCompartment::New();
    else if (compartmentType == "StationaryWater")
        additionalCompartment = anima::StationaryWaterCompartment::New();
    else if (compartmentType == "IRWater")
        additionalCompartment = anima::IsotropicRestrictedWaterCompartment::New();
    else if (compartmentType == "Stanisz")
        additionalCompartment = anima::StaniszCompartment::New();
    else if (compartmentType == "Stick")
        additionalCompartment = anima::StickCompartment::New();
    else if (compartmentType == "Zeppelin")
        additionalCompartment = anima::ZeppelinCompartment::New();
    else if (compartmentType == "Tensor")
        additionalCompartment = anima::TensorCompartment::New();
    else if (compartmentType == "NODDI")
        additionalCompartment = anima::NODDICompartment::New();
	else if (compartmentType == "DDI")
        additionalCompartment = anima::DDICompartment::New();
    else
    {
        std::string error("Unsupported compartment type: ");
        error += compartmentType;
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    return additionalCompartment;
}

} // end namespace anima
