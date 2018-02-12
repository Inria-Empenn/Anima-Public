/**
 * @file Results.cpp
 * @brief Implementation of CResults class. The class to format and saves results.
 * @author Florent Leray
 * @date 13/04/2016
 * @version 2.0
 */
#include "Results.h"

#include <string.h>
#include <stdio.h>
#include <cmath>
#include <limits>


char const*const CResults::m_ppchMeasureNameTable[eMesureLast] =
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
    "F1_score"
};


/**
   @brief    It is The default constructor.
   @details  It is private because it must be never used.
*/
CResults::CResults()
{
    for (int i=0; i<eMesureLast; ++i)
    {
        m_fResTab[i] = std::numeric_limits<float>::quiet_NaN();
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
CResults::CResults(char *pi_pchBaseFileName)
{
    for (int i=0; i<eMesureLast; ++i)
    {
        m_fResTab[i] = -1;
        m_bResActiveTab[i] = false;
    }
    if (pi_pchBaseFileName && pi_pchBaseFileName[0])
    {
        m_pchBaseOutputFileName = new char[strlen(pi_pchBaseFileName)+4+1];
        strncpy(m_pchBaseOutputFileName, pi_pchBaseFileName, strlen(pi_pchBaseFileName)+4+1);
    }
    else
    {
        m_pchBaseOutputFileName = new char[4+1];
        m_pchBaseOutputFileName[0] = 0;
    }
    m_bTxt = true;
    m_bXml = false;
    m_bScreen = false;
}

/**
   @brief    Destructor.
   @details  It prints the result on screen if necessary.
*/
CResults::~CResults()
{
    MDEL(m_pchBaseOutputFileName);
    if(m_bScreen)
    {
        for (int i=0; i<eMesureLast; ++i)
        {
            printf("%s", m_ppchMeasureNameTable[i]);
            for (int j=0; j<(20-strlen(m_ppchMeasureNameTable[i])); ++j)
            {
                printf(" ");
            }
        }
        printf("\n");
        for (int i=0; i<eMesureLast; ++i)
        {
            int iRest = 0;
            if (m_bResActiveTab[i])
            {
                iRest = printf("%f", m_fResTab[i]);
            }
            for (int j=0; j<(20-iRest); ++j)
            {
                printf(" ");
            }
        }
        printf("\n");
    }
}

/**
   @brief    It saves results on text file or xml file in function of class default settings.
   @return   True if file(s) recording are success.
*/
bool CResults::save()
{
    bool bRes = true;

    FILE *fOut = NULL;

    if (m_bTxt)
    {
        size_t len = strlen(m_pchBaseOutputFileName);
        sprintf(m_pchBaseOutputFileName+len, ".txt");
        fOut = fopen(m_pchBaseOutputFileName, "wb");
        m_pchBaseOutputFileName[len] = 0;

        if (fOut)
        {
            for (int i=0; i<eMesureLast; ++i)
            {
                if (m_bResActiveTab[i])
                {
                    bRes &= fprintf(fOut, "%f;\t", m_fResTab[i])>0;
                }
                else
                {
                    bRes &= fprintf(fOut, ";\t")>0;
                }
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
        size_t len = strlen(m_pchBaseOutputFileName);
        sprintf(m_pchBaseOutputFileName+len, ".xml");
        fOut = fopen(m_pchBaseOutputFileName, "wb");
        m_pchBaseOutputFileName[len] = 0;

        if (fOut)
        {
            char *pchImgName = strrchr(m_pchBaseOutputFileName, '/');
            if (!pchImgName)
            {
                pchImgName = strrchr(m_pchBaseOutputFileName, '\\');
                if (pchImgName)
                {
                    pchImgName++;
                }
            }
            else
            {
                pchImgName++;
            }
            if(!pchImgName)
            {
                pchImgName = m_pchBaseOutputFileName;
            }

            fprintf(fOut, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n");
            fprintf(fOut, "<image name=\"%s\">\r\n", pchImgName);
            for (int i=0; i<eMesureLast; ++i)
            {
                if (m_bResActiveTab[i])
                {
                    fprintf(fOut, "\t<measure name=\"%s\">%f</measure>\r\n", m_ppchMeasureNameTable[i], m_fResTab[i]);
                }
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
bool CResults::activeMeasurementOutput(eMesureName pi_eVal)
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
char const*const*const CResults::getMeasureNameTable()
{
    return m_ppchMeasureNameTable;
}

