#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkCommand.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <vtkPolyData.h>
#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>
#include <animaTRKReader.h>
#include <animaTRKHeaderStructure.h>

#include <iostream> // Inclure pour std::cout/cerr
#include <string>   // Inclure pour std::string
#include <exception>// Inclure pour std::exception
#include <cstdlib>  // Inclure pour EXIT_SUCCESS/FAILURE


int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("animaReadTRKInfo: Displays basic information from a .trk file.\n", ' ', "1.0");

    TCLAP::ValueArg<std::string> inArg("i", "input", "Input .trk file", true, "", "input .trk file", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        anima::TRKReader reader;
        reader.SetFileName(inArg.getValue());
        reader.Update(); // Lire le fichier et remplir m_Header et m_OutputData

        // Récupérer l'en-tête via la nouvelle méthode
        const anima::TRKHeaderStructure &header = reader.GetHeader();

        std::cout << "TRK File Information: " << inArg.getValue() << std::endl;
        std::cout << "------------------------" << std::endl;
        std::cout << "Dimensions: " << header.dim[0] << " x " << header.dim[1] << " x " << header.dim[2] << std::endl;
        std::cout << "Voxel size: " << header.voxel_size[0] << " x " << header.voxel_size[1] << " x " << header.voxel_size[2] << std::endl;
        std::cout << "Origin: " << header.origin[0] << " x " << header.origin[1] << " x " << header.origin[2] << std::endl; // Note: L'en-tête dit que ce champ n'est pas utilisé par TrackVis
        std::cout << "Number of scalars per point: " << header.n_scalars << std::endl;
        // Afficher les noms des scalaires s'il y en a
        if (header.n_scalars > 0)
        {
            std::cout << "  Scalar names: ";
            for(short i = 0; i < header.n_scalars; ++i) {
                std::cout << header.scalar_name[i] << (i == header.n_scalars - 1 ? "" : ", ");
            }
            std::cout << std::endl;
        }
        std::cout << "Number of properties per streamline: " << header.n_properties << std::endl;
         // Afficher les noms des propriétés s'il y en a
        if (header.n_properties > 0)
        {
            std::cout << "  Property names: ";
            for(short i = 0; i < header.n_properties; ++i) {
                std::cout << header.property_name[i] << (i == header.n_properties - 1 ? "" : ", ");
            }
            std::cout << std::endl;
        }
        std::cout << "Number of streamlines: " << header.n_count << (header.n_count == 0 ? " (Not stored or zero)" : "") << std::endl;
        std::cout << "Voxel order: " << header.voxel_order[0] << header.voxel_order[1] << header.voxel_order[2] << header.voxel_order[3] << std::endl;
        std::cout << "Version: " << header.version << std::endl;
        std::cout << "Header size: " << header.hdr_size << std::endl;


        // Optionnel: Afficher le polydata lu (nombre de points/cellules)
        vtkPolyData* polyData = reader.GetOutputData();
        if (polyData)
        {
             std::cout << "------------------------" << std::endl;
             std::cout << "PolyData Information:" << std::endl;
             std::cout << "Number of points: " << polyData->GetNumberOfPoints() << std::endl;
             std::cout << "Number of cells (streamlines): " << polyData->GetNumberOfCells() << std::endl;
        }


        return EXIT_SUCCESS;
    }

    catch (const itk::ExceptionObject & e)
    {
        std::cerr << "ITK Exception: " << e << std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Standard Exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "Unknown exception caught" << std::endl;
        return EXIT_FAILURE;
    }
}