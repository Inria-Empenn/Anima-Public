#include <animaTRKWriter.h>
#include <animaTRKHeaderStructure.h>
#include <fstream>

#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

namespace anima
{

void TRKWriter::Update()
{
    anima::TRKHeaderStructure headerStr;
    strcpy(headerStr.id_string,"TRACK");

    for (unsigned int i = 0;i < 3;++i)
    {
        headerStr.dim[i] = m_ReferenceImage->GetLargestPossibleRegion().GetSize(i);
        headerStr.voxel_size[i] = m_ReferenceImage->GetSpacing()[i];
        headerStr.origin[i] = 0;
    }

    // handle scalar fields
    for (unsigned int i = 0;i < 10;++i)
    {
        headerStr.scalar_name[i][0] = '\0';
        headerStr.property_name[i][0] = '\0';
    }

    unsigned int numScalarFields = m_TrackData->GetPointData()->GetNumberOfArrays();
    unsigned int numPropertyFields = m_TrackData->GetCellData()->GetNumberOfArrays();

    unsigned int pos = 0;
    std::vector <vtkDoubleArray *> scalarArrays;
    for (unsigned int i = 0;i < numScalarFields;++i)
    {
        unsigned int numComponents = m_TrackData->GetPointData()->GetArray(i)->GetNumberOfComponents();
        if (numComponents != 1)
            continue;

        std::string arrayName = m_TrackData->GetPointData()->GetArrayName(i);
        if (arrayName.size() >= 20)
            arrayName.erase(arrayName.begin() + 19);

        strcpy(headerStr.scalar_name[pos],arrayName.c_str());
        ++pos;
        scalarArrays.push_back(dynamic_cast <vtkDoubleArray *> (m_TrackData->GetPointData()->GetArray(i)));

        if (pos == 10)
            break;
    }

    headerStr.n_scalars = pos;

    pos = 0;
    std::vector <vtkDoubleArray *> propertyArrays;
    for (unsigned int i = 0;i < numPropertyFields;++i)
    {
        unsigned int numComponents = m_TrackData->GetCellData()->GetArray(i)->GetNumberOfComponents();
        if (numComponents != 1)
            continue;

        std::string arrayName = m_TrackData->GetCellData()->GetArrayName(i);
        if (arrayName.size() >= 20)
            arrayName.erase(arrayName.begin() + 19);

        strcpy(headerStr.property_name[pos],arrayName.c_str());
        propertyArrays.push_back(dynamic_cast <vtkDoubleArray *> (m_TrackData->GetCellData()->GetArray(i)));
        ++pos;
        if (pos == 10)
            break;
    }

    headerStr.n_properties = pos;

    ImageType::DirectionType direction = m_ReferenceImage->GetDirection();
    ImageType::PointType origin = m_ReferenceImage->GetOrigin();
    for (unsigned int i = 0;i < 3;++i)
    {
        for (unsigned int j = 0;j < 3;++j)
            headerStr.vox_to_ras[i][j] = direction(i,j);

        headerStr.vox_to_ras[i][3] = 0;
        headerStr.vox_to_ras[3][i] = origin[i];
    }

    headerStr.vox_to_ras[3][3] = 1;

    strcpy(headerStr.voxel_order,"RAS");
    headerStr.reserved[0] = '\0';
    headerStr.pad2[0] = '\0';
    headerStr.pad1[0] = '\0';

    // Unused flags a priori
    headerStr.invert_x = '0';
    headerStr.invert_y = '0';
    headerStr.invert_z = '0';
    headerStr.swap_xy = '0';
    headerStr.swap_yz = '0';
    headerStr.swap_zx = '0';
    for (unsigned int i = 0;i < 6;++i)
        headerStr.image_orientation_patient[i] = 0;

    headerStr.n_count = m_TrackData->GetNumberOfCells();
    headerStr.version = 2;
    headerStr.hdr_size = 1000;

    std::ofstream outFile(m_FileName,std::ios::binary);
    outFile.write((char *) &headerStr, sizeof(anima::TRKHeaderStructure));

    // handle individual tracks now : transform and store
    std::vector <float> cellData;
    std::vector <float> cellScalars(numPropertyFields);
    vnl_matrix <double> directionMatrix(4,4);

    for (unsigned int i = 0;i < 4;++i)
    {
        for (unsigned int j = 0;j < 4;++j)
            directionMatrix(i,j) = headerStr.vox_to_ras[i][j];
    }

    directionMatrix = vnl_matrix_inverse <double> (directionMatrix).as_matrix();

    for (unsigned int i = 0;i < headerStr.n_count;++i)
    {
        int cellSize = m_TrackData->GetCell(i)->GetNumberOfPoints();
        outFile.write((char *) &cellSize, sizeof(int));

        cellData.resize(3 + headerStr.n_scalars);
        for (unsigned int j = 0;j < cellSize;++j)
        {
            vtkIdType idPoint = m_TrackData->GetCell(i)->GetPointId(j);
            double *point = m_TrackData->GetPoint(m_TrackData->GetCell(i)->GetPointId(j));
            for (unsigned int k = 0;k < 3;++k)
            {
                cellData[k] = directionMatrix(k,3);
                for (unsigned int l = 0;l < 3;++l)
                    cellData[k] += directionMatrix(k,l) * point[l];
            }

                for (unsigned int k = 0;k < numScalarFields;++k)
                    cellData[k + 3] = scalarArrays[k]->GetValue(idPoint);

                outFile.write((char *) cellData.data(), (numScalarFields + 3) * sizeof(float));
        }

        if (numPropertyFields > 0)
        {
            for (unsigned int j = 0;j < numPropertyFields;++j)
                cellScalars[j] = propertyArrays[j]->GetValue(i);

            outFile.write((char *) cellScalars.data(), numPropertyFields * sizeof(float));
        }
    }

    outFile.close();
}

} // end namespace anime
