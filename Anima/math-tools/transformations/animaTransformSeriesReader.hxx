#pragma once
#include "animaTransformSeriesReader.h"

#include <tinyxml2.h>

#include <itkTransformFileReader.h>
#include <itkImageFileReader.h>

#include <itkMatrixOffsetTransformBase.h>
#include <itkStationaryVelocityFieldTransform.h>
#include <rpiDisplacementFieldTransform.h>
#include <animaVelocityUtils.h>

namespace anima
{

template <class TScalarType, unsigned int NDimensions>
TransformSeriesReader<TScalarType,NDimensions>
::TransformSeriesReader()
{
    m_OutputTransform = NULL;
    m_InvertTransform = false;
}

template <class TScalarType, unsigned int NDimensions>
TransformSeriesReader<TScalarType,NDimensions>
::~TransformSeriesReader()
{

}

template <class TScalarType, unsigned int NDimensions>
void
TransformSeriesReader<TScalarType,NDimensions>
::Update()
{
    m_OutputTransform = OutputTransformType::New();
    std::vector <TransformInformation> transformationList;

    tinyxml2::XMLDocument doc;
    tinyxml2::XMLError loadOk = doc.LoadFile(m_Input.c_str());

    if (loadOk != tinyxml2::XML_SUCCESS)
    {
        std::string error("Unable to read input summary file: ");
        error += m_Input;
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    // Read XML file into transform information vector
    tinyxml2::XMLElement *rootNode = doc.FirstChildElement( "TransformationList" );

    if (rootNode)
    {
        tinyxml2::XMLElement *trsfNode = rootNode->FirstChildElement("Transformation");

        while (trsfNode)
        {
            TransformInformation infoTrsf;
            tinyxml2::XMLElement *typeNode = trsfNode->FirstChildElement("Type");
            if (!typeNode)
                throw itk::ExceptionObject(__FILE__, __LINE__,"Type of transformation not found for a transform in the list",ITK_LOCATION);

            std::string typeStr = typeNode->GetText();
            typeStr = typeStr.substr(typeStr.find_first_not_of(" \n\r\t"));
            typeStr.erase(typeStr.find_last_not_of(" \n\r\t")+1);

            if (typeStr == "linear")
                infoTrsf.trType = LINEAR;
            else if (typeStr == "svf")
                infoTrsf.trType = SVF_FIELD;
            else if (typeStr == "dense")
                infoTrsf.trType = DENSE_FIELD;
            else
            {
                std::string error("Transformation type ");
                error += typeStr;
                error += " not supported...";
                throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
            }

            tinyxml2::XMLElement *invertNode = trsfNode->FirstChildElement("Inversion");
            if (invertNode)
            {
                std::string invertValue = invertNode->GetText();
                infoTrsf.invert = (invertValue != "0");
            }
            else
                infoTrsf.invert = false;

            if (m_InvertTransform)
                infoTrsf.invert = !infoTrsf.invert;

            tinyxml2::XMLElement *fileNode = trsfNode->FirstChildElement("Path");
            if (!fileNode)
                throw itk::ExceptionObject(__FILE__, __LINE__,"File name not found for a transform in the list",ITK_LOCATION);

            infoTrsf.fileName = fileNode->GetText();
            infoTrsf.fileName = infoTrsf.fileName.substr(infoTrsf.fileName.find_first_not_of(" \n\r\t"));
            infoTrsf.fileName.erase(infoTrsf.fileName.find_last_not_of(" \n\r\t")+1);

            transformationList.push_back(infoTrsf);

            trsfNode = trsfNode->NextSiblingElement("Transformation");
        }
    }

    // Now really load and handle global and local transform serie inversion
    // The fact that you have to apply transforms in the reverse order than the one in text file
    // is handled by the general transform

    if (m_InvertTransform)
    {
        for (int i = transformationList.size() - 1;i >= 0;--i)
        {
            switch (transformationList[i].trType)
            {
                case LINEAR:
                    this->addLinearTransformation(transformationList[i].fileName,transformationList[i].invert);
                    break;

                case SVF_FIELD:
                    this->addSVFTransformation(transformationList[i].fileName,transformationList[i].invert);
                    break;

                case DENSE_FIELD:
                default:
                    this->addDenseTransformation(transformationList[i].fileName,transformationList[i].invert);
                    break;
            }
        }
    }
    else
    {
        for (unsigned int i = 0;i < transformationList.size();++i)
        {
            switch (transformationList[i].trType)
            {
                case LINEAR:
                    this->addLinearTransformation(transformationList[i].fileName,transformationList[i].invert);
                    break;

                case SVF_FIELD:
                    this->addSVFTransformation(transformationList[i].fileName,transformationList[i].invert);
                    break;

                case DENSE_FIELD:
                default:
                    this->addDenseTransformation(transformationList[i].fileName,transformationList[i].invert);
                    break;
            }
        }
    }

    std::cout << "Loaded " << m_OutputTransform->GetNumberOfTransformsInStack() << " transformations from transform list file: " << m_Input << std::endl;
}

template <class TScalarType, unsigned int NDimensions>
void
TransformSeriesReader<TScalarType,NDimensions>
::addLinearTransformation(std::string &fileName,bool invert)
{
    typedef itk::MatrixOffsetTransformBase <TScalarType,NDimensions> MatrixTransformType;
    typedef typename MatrixTransformType::Pointer MatrixTransformPointer;

    itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
    reader->SetFileName(fileName);
    reader->Update();

    const itk::TransformFileReader::TransformListType *trsfList = reader->GetTransformList();
    itk::TransformFileReader::TransformListType::const_iterator tr_it = trsfList->begin();

    MatrixTransformPointer trsf = dynamic_cast <MatrixTransformType *> ((*tr_it).GetPointer());

    if (invert)
    {
        MatrixTransformPointer tmpInvert = MatrixTransformType::New();
        trsf->GetInverse(tmpInvert);

        trsf = tmpInvert;
    }

    m_OutputTransform->InsertTransform(trsf.GetPointer());
}

template <class TScalarType, unsigned int NDimensions>
void
TransformSeriesReader<TScalarType,NDimensions>
::addDenseTransformation(std::string &fileName,bool invert)
{
    typedef rpi::DisplacementFieldTransform <TScalarType,NDimensions> DenseTransformType;
    typedef typename DenseTransformType::Pointer DenseTransformPointer;
    typedef typename DenseTransformType::VectorFieldType DisplacementFieldType;

    typedef itk::ImageFileReader<DisplacementFieldType> DispReaderType;
    typename DispReaderType::Pointer trReader = DispReaderType::New();
    trReader->SetFileName(fileName);
    trReader->Update();

    DenseTransformPointer dispTrsf = DenseTransformType::New();
    dispTrsf->SetParametersAsVectorField(trReader->GetOutput());

    if (invert)
    {
        DenseTransformPointer tmpInvert = DenseTransformType::New();
        dispTrsf->GetInverse(tmpInvert);

        dispTrsf = tmpInvert;
    }

    m_OutputTransform->InsertTransform(dispTrsf.GetPointer());
}

template <class TScalarType, unsigned int NDimensions>
void
TransformSeriesReader<TScalarType,NDimensions>
::addSVFTransformation(std::string &fileName,bool invert)
{
    typedef itk::StationaryVelocityFieldTransform <TScalarType,NDimensions> SVFTransformType;
    typedef typename SVFTransformType::Pointer SVFTransformPointer;
    typedef typename SVFTransformType::VectorFieldType VelocityFieldType;

    typedef rpi::DisplacementFieldTransform <TScalarType,NDimensions> DenseTransformType;
    typedef typename DenseTransformType::Pointer DenseTransformPointer;

    typedef itk::ImageFileReader<VelocityFieldType> SVFReaderType;
    typename SVFReaderType::Pointer trReader = SVFReaderType::New();
    trReader->SetFileName(fileName);
    trReader->Update();

    SVFTransformPointer svfPointer = SVFTransformType::New();
    svfPointer->SetParametersAsVectorField(trReader->GetOutput());

    DenseTransformPointer dispTrsf = DenseTransformType::New();
    anima::GetSVFExponential(svfPointer.GetPointer(),dispTrsf.GetPointer(),invert);

    m_OutputTransform->InsertTransform(dispTrsf.GetPointer());
}

} // end of namespace anima
