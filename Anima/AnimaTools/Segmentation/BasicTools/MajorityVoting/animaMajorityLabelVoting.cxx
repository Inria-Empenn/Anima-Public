#include <animaMajorityLabelVotingImageFilter.h>
#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i","input-images","Input images",true,"","Label input images (list or 4D image)",cmd);
    TCLAP::ValueArg<std::string> consensusImageArg("o","consensus-image","consensus label image",true,"","consensus image",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <unsigned int, 3> ImageType;
    typedef anima::MajorityLabelVotingImageFilter <unsigned int> VotingFilterType;

    VotingFilterType::Pointer votingFilter = VotingFilterType::New();
    anima::setMultipleImageFilterInputsFromFileName <ImageType,VotingFilterType> (inputArg.getValue(),votingFilter);

    votingFilter->Update();

    anima::writeImage <ImageType> (consensusImageArg.getValue(), votingFilter->GetOutput());

    return EXIT_SUCCESS;
}
