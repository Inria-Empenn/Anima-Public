#include <tclap/CmdLine.h>

#include <string>
#include <fstream>

#include <itkTransformFileReader.h>
#include <itkImageFileReader.h>

int main(int ac, const char **av)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::MultiArg<std::string> inputArg("i", "input", "multiple input transform names",true,"input files",cmd);
    TCLAP::MultiArg<unsigned int> invertArg("I", "invert", "multiple invert transform arguments",false,"invert inputs",cmd);
    
    TCLAP::ValueArg<std::string> outputArg("o","output","output xml filename",true,"","output filename",cmd);
    
    TCLAP::SwitchArg denseArg("D","dense","Non linear transfroms are dense fields (default: they are SVFs)",cmd,false);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    std::vector <std::string> inputNames = inputArg.getValue();
    std::vector <unsigned int> inputInverts = invertArg.getValue();
    
    bool useInvert = (inputInverts.size() == inputNames.size());
    
    std::ofstream outputFile(outputArg.getValue().c_str());
    outputFile << "<?xml version=\"1.0\"?>" << std::endl;
    outputFile << "<TransformationList>" << std::endl;
    
    for (unsigned int i = 0;i < inputNames.size();++i)
    {
        // Try first to read as a linear transform
        bool linearSuccess = true;
        try
        {
            itk::TransformFileReader::Pointer trReader = itk::TransformFileReader::New();
            trReader->SetFileName(inputNames[i]);
            
            trReader->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            linearSuccess = false;
        }
        
        if (linearSuccess)
        {
            outputFile << "<Transformation>" << std::endl;
            outputFile << "<Type>linear</Type>" << std::endl;
            
            if (useInvert)
                outputFile << "<Inversion>" << inputInverts[i] << "</Inversion>" << std::endl;
            else
                outputFile << "<Inversion>0</Inversion>" << std::endl;
            
            outputFile << "<Path>" << inputNames[i] << "</Path>" << std::endl;
            outputFile << "</Transformation>" << std::endl;
            
            continue;
        }
        
        // Then try non linear transform
        bool nonLinearSuccess = true;
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputNames[i].c_str(), itk::ImageIOFactory::ReadMode);

        if (!imageIO)
            nonLinearSuccess = false;
        
        if (nonLinearSuccess)
        {
            outputFile << "<Transformation>" << std::endl;
            
            if (denseArg.isSet())
                outputFile << "<Type>dense</Type>" << std::endl;
            else
                outputFile << "<Type>svf</Type>" << std::endl;
            
            if (useInvert)
                outputFile << "<Inversion>" << inputInverts[i] << "</Inversion>" << std::endl;
            else
                outputFile << "<Inversion>0</Inversion>" << std::endl;
            
            outputFile << "<Path>" << inputNames[i] << "</Path>" << std::endl;
            outputFile << "</Transformation>" << std::endl;
            
            continue;
        }
        
        // else error and quit
        std::cerr << inputNames[i] << " is not a readable transform, exiting" << std::endl;
        return 1;
    }
    
    outputFile << "</TransformationList>" << std::endl;
    outputFile.close();
    
    std::cout << "Generated an XML file " << outputArg.getValue() << " with " << inputNames.size() << " transformations..." << std::endl;

    return 0;
}
