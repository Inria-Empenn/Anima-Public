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
    double m_fSensitivity;
    double m_fSpecificity;
    double m_fPPV;
    double m_fNPV;
    double m_fDice;
    double m_fJaccard;
    double m_fRVE;

    //////////////////////////////////////////////////////////////////////////
    // Distances metrics
    double m_fHausdorffDist;
    double m_fMeanDist;
    double m_fAverageDist;

    //////////////////////////////////////////////////////////////////////////
    // Detection scores
    double m_fPPVL;
    double m_fSensL;
    double m_fF1;

    //////////////////////////////////////////////////////////////////////////
    // general informations
    int m_iNbThreads;     /*<! number of thread used by processing. */
    std::string m_pchOutBase;   /*<! base name used for output file. */

    //////////////////////////////////////////////////////////////////////////
    // Detection scores parameters
    double m_fDetectionLesionMinVolume;
    double m_fTPLMinOverlapRatio;
    double m_fTPLMaxFalsePositiveRatio;
    double m_fTPLMaxFalsePositiveRatioModerator;

    std::string m_oStrInImage;   /*<! Path of Image to test. */
    std::string m_oStrRefImage;  /*<! Path of reference Image. */
    std::string m_oStrBaseOut;   /*<! Base name for output results file. */
};

} // end namespace anima
