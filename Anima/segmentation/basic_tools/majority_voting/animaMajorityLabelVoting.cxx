#include <itkImageRegionConstIterator.h>
#include <itkImageLinearConstIteratorWithIndex.h>
#include <itkImage.h>

#include <animaReadWriteFunctions.h>

#include <limits.h>
#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>

#include <map>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i","inputimage","Input image",true,"","Input image (4D image, concatenation of multiple 3D label images)",cmd);
    TCLAP::ValueArg<std::string> consensusImageArg("o","consensusImage","consensus label image",false,"","consensus image",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <int,4> Image4DType;
    typedef itk::Image <int,3> Image3DType;

    typedef itk::ImageLinearConstIteratorWithIndex <Image4DType> IteratorType;

    Image3DType::IndexType index3D;
    Image4DType::IndexType index4D;

    Image4DType::Pointer inputImage = anima::readImage <Image4DType> (inputArg.getValue());

    Image3DType::PointType outputOrigin;
    Image3DType::SpacingType outputSpacing;
    Image3DType::DirectionType outputDirection;
    Image3DType::SizeType outputSize;
    Image3DType::RegionType outputRegion;

    for (unsigned int d = 0; d < 3; ++d)
    {
        outputOrigin[d] = inputImage->GetOrigin()[d]; 
        outputSpacing[d] = inputImage->GetSpacing()[d]; 
        outputSize[d] = inputImage->GetLargestPossibleRegion().GetSize()[d];
        for(unsigned int i = 0; i < 3; ++i)
            outputDirection[d][i] = inputImage->GetDirection()[d][i];     
    } 
    outputRegion.SetSize(outputSize);

    Image3DType::Pointer consensusImage = Image3DType::New();
    if(consensusImageArg.getValue() != "")
    {
        consensusImage->Initialize();
        consensusImage->SetOrigin(outputOrigin);
        consensusImage->SetSpacing(outputSpacing);
        consensusImage->SetDirection(outputDirection);
        consensusImage->SetRegions(outputRegion);
        consensusImage->Allocate();
    }

    Image4DType::RegionType tmpRegionInputImage = inputImage->GetLargestPossibleRegion();
    const unsigned int timeLength = tmpRegionInputImage.GetSize()[3];

    IteratorType it(inputImage,tmpRegionInputImage);
    it.SetDirection(3);
    it.GoToBegin();


    std::vector <int> tmpVec (timeLength);

    while( !it.IsAtEnd() )
    {
        it.GoToBeginOfLine();
        index4D = it.GetIndex();

        for (unsigned int i = 0;i < timeLength;++i)
        {
            tmpVec[i] = it.Get();
            ++it;
        }

        int consensusLabel = mostFrequent(tmpVec, timeLength);

        index3D[0] = index4D[0];
        index3D[1] = index4D[1];
        index3D[2] = index4D[2];

        if(consensusImageArg.getValue() != "")
            consensusImage->SetPixel( index3D, consensusLabel );

        it.NextLine();
    }
    
    if(consensusImageArg.getValue() != "")
        anima::writeImage <Image3DType> (consensusImageArg.getValue(), consensusImage);

    return EXIT_SUCCESS;
}


int mostFrequent(int arr[], unsigned int n) 
{ 
    sort(arr, arr + n); 
  
    int max_count = 1, res = arr[0], curr_count = 1; 
    for (int i = 1; i < n; i++) { 
        if (arr[i] == arr[i - 1]) 
            curr_count++; 
        else { 
            if (curr_count > max_count) { 
                max_count = curr_count; 
                res = arr[i - 1]; 
            } 
            curr_count = 1; 
        } 
    } 
  
    if (curr_count > max_count) 
    { 
        max_count = curr_count; 
        res = arr[n - 1]; 
    } 
  
    return res; 
} 

