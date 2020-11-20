#include <animaMCMFileReader.h>
#include <tinyxml2.h>
#include <itkObjectFactoryBase.h>

namespace anima
{

itk::IOComponentEnum GetMCMComponentType(std::string fileName)
{
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLError loadOk = doc.LoadFile(fileName.c_str());

    std::replace(fileName.begin(),fileName.end(),'\\','/');
    std::string basePath, baseName;
    std::size_t lastSlashPos = fileName.find_last_of("/");

    if (lastSlashPos != std::string::npos)
        basePath.append(fileName.begin(),fileName.begin() + lastSlashPos);

    if (basePath.length() != 0)
        basePath = basePath + "/";

    std::size_t lastDotPos = fileName.find_last_of(".");
    baseName.append(fileName.begin() + lastSlashPos, fileName.begin() + lastDotPos);

    if (loadOk != tinyxml2::XML_SUCCESS)
        return itk::IOComponentEnum::UNKNOWNCOMPONENTTYPE;

    // Read XML file into transform information vector
    tinyxml2::XMLElement *modelNode = doc.FirstChildElement( "Model" );
    if (!modelNode)
        return itk::IOComponentEnum::UNKNOWNCOMPONENTTYPE;

    tinyxml2::XMLElement *weightsNode = modelNode->FirstChildElement( "Weights" );

    if (!weightsNode)
        return itk::IOComponentEnum::UNKNOWNCOMPONENTTYPE;

    std::string weightsFileName = basePath + baseName + "/";
    weightsFileName += weightsNode->GetText();

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(weightsFileName.c_str(),
                                                                           itk::IOFileModeEnum::ReadMode);

    if (!imageIO)
        return itk::IOComponentEnum::UNKNOWNCOMPONENTTYPE;

    imageIO->SetFileName(weightsFileName.c_str());

    try
    {
        imageIO->ReadImageInformation();
    }
    catch(itk::ExceptionObject &e)
    {
        return itk::IOComponentEnum::UNKNOWNCOMPONENTTYPE;
    }

    return imageIO->GetComponentType();
}

} // end namespace anima
