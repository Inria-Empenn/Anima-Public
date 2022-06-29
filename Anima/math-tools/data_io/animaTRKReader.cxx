#include <animaTRKReader.h>
#include <animaTRKHeaderStructure.h>
#include <vnl_matrix.h>
#include <fstream>

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include <itkMacro.h>

namespace anima
{

void TRKReader::Update()
{
    anima::TRKHeaderStructure headerStr;

    std::ifstream inFile(m_FileName,std::ios::binary);
    if (!inFile.is_open())
        throw itk::ExceptionObject(__FILE__, __LINE__,"Unable to open file " + m_FileName,ITK_LOCATION);

    inFile.read((char *) &headerStr, sizeof(anima::TRKHeaderStructure));

    if (headerStr.version != 2)
        throw itk::ExceptionObject(__FILE__, __LINE__,"TRK reader only supports version 2",ITK_LOCATION);

    // Now create polydata and fill it
    m_OutputData = vtkSmartPointer <vtkPolyData>::New();
    m_OutputData->Initialize();
    m_OutputData->Allocate();

    vtkSmartPointer <vtkPoints> myPoints = vtkSmartPointer <vtkPoints>::New();
    std::vector < vtkSmartPointer <vtkDoubleArray> > scalarArrays(headerStr.n_scalars);
    for (unsigned int i = 0;i < headerStr.n_scalars;++i)
    {
        scalarArrays[i] = vtkSmartPointer <vtkDoubleArray>::New();
        scalarArrays[i]->SetNumberOfComponents(1);
        scalarArrays[i]->SetName(headerStr.scalar_name[i]);
    }

    std::vector < vtkSmartPointer <vtkDoubleArray> > cellArrays(headerStr.n_properties);
    for (unsigned int i = 0;i < headerStr.n_properties;++i)
    {
        cellArrays[i] = vtkSmartPointer <vtkDoubleArray>::New();
        cellArrays[i]->SetNumberOfComponents(1);
        cellArrays[i]->SetName(headerStr.property_name[i]);
    }

    unsigned int nCells = headerStr.n_count;
    std::vector <float> pointValues(3 + headerStr.n_scalars);
    std::vector <float> cellScalars(headerStr.n_properties);
    vnl_matrix <double> vox_to_ras(4,4);
    for (unsigned int i = 0;i < 4;++i)
    {
        for (unsigned int j = 0;j < 4;++j)
            vox_to_ras(i,j) = headerStr.vox_to_ras[i][j];
    }

    for (unsigned int i = 0;i < nCells;++i)
    {
        int npts;
        inFile.read((char *) &npts, sizeof(int));

        vtkIdType* ids = new vtkIdType[npts];

        for (unsigned int j = 0;j < npts;++j)
        {
            inFile.read((char *) pointValues.data(), (3 + headerStr.n_scalars) * sizeof(float));

            double xValue = vox_to_ras(0,3);
            double yValue = vox_to_ras(1,3);
            double zValue = vox_to_ras(2,3);

            for (unsigned int k = 0;k < 3;++k)
            {
                xValue += vox_to_ras(0,k) * pointValues[k];
                yValue += vox_to_ras(1,k) * pointValues[k];
                zValue += vox_to_ras(2,k) * pointValues[k];
            }

            ids[j] = myPoints->InsertNextPoint(xValue, yValue, zValue);
            for (unsigned int k = 0;k < headerStr.n_scalars;++k)
                scalarArrays[k]->InsertNextValue(pointValues[3 + k]);
        }

        inFile.read((char *) cellScalars.data(), headerStr.n_properties * sizeof(float));
        for (unsigned int k = 0;k < headerStr.n_properties;++k)
            cellArrays[k]->InsertNextValue(cellScalars[k]);

        m_OutputData->InsertNextCell (VTK_POLY_LINE, npts, ids);
        delete[] ids;
    }

    inFile.close();

    m_OutputData->SetPoints(myPoints);
    for (unsigned int k = 0;k < headerStr.n_scalars;++k)
        m_OutputData->GetPointData()->AddArray(scalarArrays[k]);
    for (unsigned int k = 0;k < headerStr.n_properties;++k)
        m_OutputData->GetCellData()->AddArray(cellArrays[k]);
}

} // end namespace anime
