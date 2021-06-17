#include "animaSegPerfResults.h"

#include <string.h>
#include <stdio.h>
#include <cmath>
#include <limits>

namespace anima
{

char const*const SegPerfResults::m_ppchMeasureNameTable[eMesureLast] =
{
    "Jaccard",
    "Dice",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NPV",
    "RelativeVolumeError",
    "HausdorffDistance",
    "ContourMeanDistance",
    "SurfaceDistance",
    "PPVL",
    "SensL",
    "F1_score",
    "NbTestedLesions",
    "VolTestedLesions"
};

/**
   @brief    It is The default constructor.
   @details  It is private because it must be never used.
*/
SegPerfResults::SegPerfResults()
{
    for (int i=0; i<eMesureLast; ++i)
    {
        m_fResTab[i] = std::numeric_limits<double>::quiet_NaN();
    }
    m_bTxt = true;
    m_bXml = false;
    m_bScreen = false;
    m_pchBaseOutputFileName = new char[4+1];
}

/**
   @brief    It is the constructor to used.
   @param    [in] pi_pchBaseFileName Name of the file to evaluate.
*/
SegPerfResults::SegPerfResults(std::string &pi_pchBaseFileName)
{
    for (int i=0; i<eMesureLast; ++i)
    {
        m_fResTab[i] = -1;
        m_bResActiveTab[i] = false;
    }

    if (pi_pchBaseFileName != "")
        m_pchBaseOutputFileName = pi_pchBaseFileName;

    m_bTxt = true;
    m_bXml = false;
    m_bScreen = false;
}

/**
   @brief    Destructor.
   @details  It prints the result on screen if necessary.
*/
SegPerfResults::~SegPerfResults()
{
    if(m_bScreen)
    {
        for (int i=0; i<eMesureLast; ++i)
        {
            printf("%s", m_ppchMeasureNameTable[i]);
            for (int j=0; j<(20-strlen(m_ppchMeasureNameTable[i])); ++j)
                printf(" ");
        }
        printf("\n");

        for (int i=0; i<eMesureLast; ++i)
        {
            int iRest = 0;
            if (m_bResActiveTab[i])
                iRest = printf("%f", m_fResTab[i]);

            for (int j=0; j<(20-iRest); ++j)
                printf(" ");
        }
        printf("\n");
    }
}

/**
   @brief    It saves results on text file or xml file in function of class default settings.
   @return   True if file(s) recording are success.
*/
bool SegPerfResults::save()
{
    bool bRes = true;

    FILE *fOut = NULL;

    if (m_bTxt)
    {
        std::string outFileName = m_pchBaseOutputFileName + ".txt";
        fOut = fopen(outFileName.c_str(), "wb");

        if (fOut)
        {
            for (int i=0; i<eMesureLast; ++i)
            {
                if (m_bResActiveTab[i])
                    bRes &= fprintf(fOut, "%f;\t", m_fResTab[i])>0;
                else
                    bRes &= fprintf(fOut, ";\t")>0;
            }
            bRes &= fprintf(fOut, "\r\n")>0;
            fclose(fOut);
        }
        else
        {
            bRes = false;
        }
    }

    if (m_bXml)
    {
        std::string outFileName = m_pchBaseOutputFileName + ".xml";
        fOut = fopen(outFileName.c_str(), "wb");

        if (fOut)
        {
            char *tmpStr = const_cast <char *> (m_pchBaseOutputFileName.c_str());
            char *pchImgName = strrchr(tmpStr, '/');
            if (!pchImgName)
            {
                pchImgName = strrchr(tmpStr, '\\');
                if (pchImgName)
                    pchImgName++;
            }
            else
                pchImgName++;

            if(!pchImgName)
                pchImgName = tmpStr;

            fprintf(fOut, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n");
            fprintf(fOut, "<image name=\"%s\">\r\n", pchImgName);
            for (int i=0; i<eMesureLast; ++i)
            {
                if (m_bResActiveTab[i])
                    fprintf(fOut, "\t<measure name=\"%s\">%f</measure>\r\n", m_ppchMeasureNameTable[i], m_fResTab[i]);
            }
            fprintf(fOut, "</image>\r\n");
            fclose(fOut);
        }
        else
        {
            bRes = false;
        }
    }

    return bRes;
}

/**
@brief    It active the saving of one specific measure. If it set twice time the effect is inverted.
@return   True if file(s) recording are success.
*/
bool SegPerfResults::activeMeasurementOutput(eMesureName pi_eVal)
{
    bool bRes = false;

    if (pi_eVal < eMesureLast)
    {
        m_bResActiveTab[pi_eVal] = !m_bResActiveTab[pi_eVal];
        bRes = m_bResActiveTab[pi_eVal];
    }

    return bRes;
}

/**
@brief    Get the list of all Measures available.
@return   A constant ordered list of measures names.
*/
char const*const*const SegPerfResults::getMeasureNameTable()
{
    return m_ppchMeasureNameTable;
}

} // end namespace anima
