// filepath: /home/ndecaux/Git/Anima/src/Anima/math-tools/data_io/animaTRKReader.h
#pragma once

#include <libAnimaCoreExport.h>
#include <animaTRKHeaderStructure.h> // Inclure la définition de la structure
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <string> // Inclure pour std::string

namespace anima
{

class LIBANIMACORE_EXPORT TRKReader
{
public:
    TRKReader()
    {
        m_FileName = "";
        // Initialiser la structure header si nécessaire, par exemple avec des zéros
        memset(&m_Header, 0, sizeof(TRKHeaderStructure));
    }

    ~TRKReader() {}

    void SetFileName(const std::string &name) {m_FileName = name;} 
    vtkPolyData *GetOutputData() {return m_OutputData;}
    const TRKHeaderStructure& GetHeader() const { return m_Header; } // Getter pour l'en-tête

    void Update();

private:
    vtkSmartPointer <vtkPolyData> m_OutputData;
    std::string m_FileName;
    TRKHeaderStructure m_Header; // Membre pour stocker l'en-tête
};

} // end namespace anima
