#include <tclap/CmdLine.h>
#include <sstream>


#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkTransform.h>



int main(int ac, const char** av)
{

    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
     
    // Scale
    TCLAP::ValueArg<double> 
            tArg("t",
                 "t_value",
                 "Factor of scaling for the z axis.",
                 false,
                 1,
                 "double",
                 cmd
                 );
    
    TCLAP::ValueArg<double> 
            sArg("s",
                 "s_value",
                 "Factor of scaling for the s axis.",
                 false,
                 1,
                 "double",
                 cmd
                  );
    TCLAP::ValueArg<double> 
            rArg("r",
                 "r_value",
                 "Factor of scaling for the r axis.",
                 false,
                 1,
                 "double",
                 cmd
                 );
    
    // Rotation
    TCLAP::ValueArg<double> 
            wArg("w",
                 "w_value",
                 "Angle of the translation aroud the z axis.",
                 false,
                 0,
                 "double",
                 cmd
                 );
    
    TCLAP::ValueArg<double> 
            vArg("v",
                 "v_value",
                 "Angle of the translation aroud the y axis",
                 false,
                 0,
                 "double",
                 cmd
                  );
    TCLAP::ValueArg<double> 
            uArg("u",
                 "u_value",
                 "Angle of the translation aroud the x axis",
                 false,
                 0,
                 "double",
                 cmd
                 );
    
    
    // Translation
    TCLAP::ValueArg<double> 
            zArg("z",
                 "z_value",
                 "Z value of the translation.",
                 false,
                 0,
                 "double",
                 cmd
                 );
    
    TCLAP::ValueArg<double> 
            yArg("y",
                 "y_value",
                 "Z value of the translation",
                 false,
                 0,
                 "double",
                 cmd
                  );
    TCLAP::ValueArg<double> 
            xArg("x",
                 "x_value",
                 "X value of the translation",
                 false,
                 0,
                 "double",
                 cmd
                 );

    
    
    // Files
    TCLAP::ValueArg<std::string> 
            outputArg("o",
                      "output",
                      "the registered moving object",
                      true,
                      "",
                      "Filename",
                      cmd
                      );
    TCLAP::ValueArg<std::string> 
            inputArg("i",
                     "reference",
                     "The reference object",
                     true,
                     "",
                     "Filename",
                     cmd
                     );
    
    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    vtkPolyData *data = vtkPolyData::New();
    vtkPoints *points = vtkPoints::New();
    vtkPoints *newPoints = vtkPoints::New();
    

    // Read from input
    reader->SetFileName(inputArg.getValue().c_str());
    reader->Update();
    data = reader->GetOutput();
    points = data->GetPoints();
    
    vtkTransform *transform = vtkTransform::New();
    // Set translation
    transform->Translate(xArg.getValue(), yArg.getValue(), zArg.getValue());

    // Set Rotation
    transform->RotateX(uArg.getValue());
    transform->RotateY(vArg.getValue());
    transform->RotateZ(wArg.getValue());
    
    // Set Scale
    transform->Scale(rArg.getValue(), sArg.getValue(), tArg.getValue());
    
    // Apply to output
    std::cout<<"transformation matrix : ";
    transform->GetMatrix()->Print(std::cout);
    transform->TransformPoints(points, newPoints);
    data->SetPoints(newPoints);
        
    // Write to output
    std::string outFileName = outputArg.getValue();
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(outFileName.c_str());
    std::cout<< "write in : " << outFileName.c_str() <<std::endl;
    writer->SetInputData(data);
    writer->Update();
}
    
