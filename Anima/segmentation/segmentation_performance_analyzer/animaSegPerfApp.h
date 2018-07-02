#pragma once

#include <animaSegPerfResults.h>
#include <animaSegPerfCAnalyzer.h>

namespace anima
{

/**
* @class SegPerfApp
* @brief Main class to structure application and handle command line options.
*/
class SegPerfApp
{
public:
    SegPerfApp(void);
    ~SegPerfApp(void);

    bool init(int argc, char *argv[]);
    bool checkParamsCoherence();
    void checkOutputCoherence();
    void prepareOutput();
    void play();

    static void about();

protected:
    void processAnalyze(SegPerfCAnalyzer &pi_oAnalyzer, int pi_iIndex);
    void storeMetricsAndMarks(SegPerfResults &pi_roRes);
    long writeStoredMetricsAndMarks(SegPerfResults &pi_roRes);

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
};

} // end namespace anima
