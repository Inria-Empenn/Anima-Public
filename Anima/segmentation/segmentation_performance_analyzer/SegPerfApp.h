#pragma once

#include <Results.h>
#include <Analyzer.h>

/**
* @class CSegPerfApp
* @brief Main class to structure application and handle command line options.
*/
class CSegPerfApp
{
private:
    //////////////////////////////////////////////////////////////////////////
    // Output way
    bool m_bTxt;
    bool m_bXml;
    bool m_bScreen;

    //////////////////////////////////////////////////////////////////////////
    // Group of metrics enable
    bool m_bSegmentationEvaluation;
    bool m_bAdvancedEvaluation;
    bool m_bSurfaceEvaluation;
    bool m_bLesionsDetectionEvaluation;

    //////////////////////////////////////////////////////////////////////////
    // Lesions specific metrics
    float m_fSensitivity;
    float m_fSpecificity;
    float m_fPPV;
    float m_fNPV;
    float m_fDice;
    float m_fJaccard;
    float m_fRVE;

    //////////////////////////////////////////////////////////////////////////
    // Distances metrics
    float m_fHausdorffDist;
    float m_fMeanDist;
    float m_fAverageDist;

    //////////////////////////////////////////////////////////////////////////
    // Detection scores
    float m_fPPVL;
    float m_fSensL;
    float m_fF1;

    //////////////////////////////////////////////////////////////////////////
    // general informations
    int m_iNbThreads;     /*<! number of thread used by processing. */
    std::string m_pchOutBase;   /*<! base name used for output file. */

    //////////////////////////////////////////////////////////////////////////
    // Detection scores parameters
    float m_fDetectionLesionMinVolume;
    float m_fTPLMinOverlapRatio;
    float m_fTPLMaxFalsePositiveRatio;
    float m_fTPLMaxFalsePositiveRatioModerator;

    std::string m_oStrInImage;   /*<! Path of Image to test. */
    std::string m_oStrRefImage;  /*<! Path of reference Image. */
    std::string m_oStrBaseOut;   /*<! Base name for output results file. */

public:
    CSegPerfApp(void);
    ~CSegPerfApp(void);

    bool init(int argc, char *argv[]);
    bool checkParamsCoherence();
    void checkOutputCoherence();
    void prepareOutput();
    void play();

    static void about();

private:
    void processAnalyze(CAnalyzer&pi_oAnalyzer, int pi_iIndex);
    void storeMetricsAndMarks(CResults&pi_roRes);
    long writeStoredMetricsAndMarks(CResults&pi_roRes);
};
