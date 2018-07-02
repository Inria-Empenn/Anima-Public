#include "SegPerfApp.h"

#include <tclap/CmdLine.h>
#include <iostream>
#include <cmath>
#include <limits>

CSegPerfApp::CSegPerfApp(void)
{
    //////////////////////////////////////////////////////////////////////////
    // Output way
    m_bTxt = true;
    m_bXml = false;
    m_bScreen = false;

    //////////////////////////////////////////////////////////////////////////
    // Group of metrics enable
    m_bSegmentationEvaluation = true;
    m_bAdvancedEvaluation = false;
    m_bSurfaceEvaluation = true;
    m_bLesionsDetectionEvaluation = true;

    //////////////////////////////////////////////////////////////////////////
    // Lesions specific metrics
    m_fSensitivity = std::numeric_limits<float>::quiet_NaN();
    m_fSpecificity = std::numeric_limits<float>::quiet_NaN();
    m_fPPV = std::numeric_limits<float>::quiet_NaN();
    m_fNPV = std::numeric_limits<float>::quiet_NaN();
    m_fDice = std::numeric_limits<float>::quiet_NaN();
    m_fJaccard = std::numeric_limits<float>::quiet_NaN();
    m_fRVE = std::numeric_limits<float>::quiet_NaN();

    //////////////////////////////////////////////////////////////////////////
    // Distances metrics
    m_fHausdorffDist = std::numeric_limits<float>::quiet_NaN();
    m_fMeanDist = std::numeric_limits<float>::quiet_NaN();
    m_fAverageDist = std::numeric_limits<float>::quiet_NaN();

    //////////////////////////////////////////////////////////////////////////
    // Detection scores
    m_fPPVL = std::numeric_limits<float>::quiet_NaN();
    m_fSensL = std::numeric_limits<float>::quiet_NaN();
    m_fF1 = std::numeric_limits<float>::quiet_NaN();

    //////////////////////////////////////////////////////////////////////////
    // general informations
    m_iNbThreads = 0;      /*<! number of thread used by processing. */
    m_pchOutBase = "";   /*<! base name used for output file. */

    //////////////////////////////////////////////////////////////////////////
    // Detection scores parameters
    m_fDetectionLesionMinVolume = 3.0;
    m_fTPLMinOverlapRatio = 0.1;
    m_fTPLMaxFalsePositiveRatio = 0.7;
    m_fTPLMaxFalsePositiveRatioModerator = 0.65;
}

CSegPerfApp::~CSegPerfApp(void)
{
}

/**
   @brief    This method set the application behavior thanks to command line arguments parsing.
   @param    [in] argc is the number of command line arguments.
   @param    [in] argv is the table  of command line arguments.
   @details  The main function must delegate to it command line arguments. this method will parse and set different option values.
*/
bool CSegPerfApp::init(int argc, char *argv[])
{
    // Define the command line object.
    TCLAP::CmdLine cmd("Tools to analyze segmentation performances by comparison", ' ', ANIMA_VERSION);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> oArgInputImg("i", "input", "Input image.", true, "", "string", cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> oArgRefImg("r", "ref", "Reference image to compare input image.", true, "", "string", cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> oArgBaseOutputName("o", "outputBase", "Base name for output files", false, "", "string", cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<int> oArgNbThreads("t", "threads", "Number of threads", false, 0, "int", cmd);

    // Define a value argument and add it to the command line.
    TCLAP::SwitchArg oArgAbout("A", "About", "Details on output metrics", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg oArgSwitchText("T", "text", "Stores results into a text file.", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg oArgSwitchXml("X", "xml", "Stores results into a xml file.", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg oArgSwitchScreen("S", "screen", "Print results on the screen when program ended.", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg oArgSegEval("s", "SegmentationEvaluationMetrics", "Compute metrics to evaluate a segmentation.", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg oArgSwitchAdvancedEvaluation("a", "advancedEvaluation", "Compute results for each cluster (intra-lesion results)", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg oArgSwitchDetectionEval("l", "LesionDetectionMetrics", "Compute metrics to evaluate the detection of lesions along to a segmentation.", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg oArgSwitchSurfaceDist("d", "SurfaceEvaluation", "Surface distances evaluation.", cmd, false);

    // Define a switch and add it to the command line.
    TCLAP::ValueArg<float> oArgSwitchDetectionLesionMinVolume("v", "MinLesionVolume", "Min volume of lesion for \"Lesions detection metrics\" in mm^3 (default 3mm^3).", false, 3.00, "float", cmd);

    // Define a switch and add it to the command line.
    TCLAP::ValueArg<float> oArgSwitchTPLMinOverlapRatio("x", "MinOverlapRatio", "Minimum overlap ratio to say if a lesion of the GT is detected. (default 0.10)", false, 0.10, "float", cmd);

    // Define a switch and add it to the command line.
    TCLAP::ValueArg<float> oArgSwitchTPLMaxFalsePositiveRatio("y", "MaxFalsePositiveRatio", "Maximum of false positive ratio to limit the detection of a lesion in GT if a lesion in the image is too big. (default 0.7)", false, 0.70, "float", cmd);

    // Define a switch and add it to the command line.
    TCLAP::ValueArg<float> oArgSwitchTPLMaxFalsePositiveRatioModerator("z", "MaxFalsePositiveRatioModerator", "Percentage of the regions overlapping the tested lesion is not too much outside of this lesion. (default 0.65)", false, 0.65, "float", cmd);

    // Parse the args.
    cmd.parse( argc, argv );

    if (oArgAbout.isSet())
    {
        about();
        return false;
    }

    m_oStrInImage  = oArgInputImg.getValue();
    m_oStrRefImage = oArgRefImg.getValue();
    m_oStrBaseOut  = oArgBaseOutputName.getValue();

    m_bTxt         = oArgSwitchText.getValue();
    m_bXml         = oArgSwitchXml.getValue();
    m_bScreen      = oArgSwitchScreen.getValue();

    m_iNbThreads   = oArgNbThreads.getValue();

    m_bSegmentationEvaluation = oArgSegEval.getValue();
    m_bAdvancedEvaluation = oArgSwitchAdvancedEvaluation.getValue();
    m_bSurfaceEvaluation = oArgSwitchSurfaceDist.getValue();
    m_bLesionsDetectionEvaluation = oArgSwitchDetectionEval.getValue();

    m_fDetectionLesionMinVolume = oArgSwitchDetectionLesionMinVolume.getValue();
    m_fTPLMinOverlapRatio = oArgSwitchTPLMinOverlapRatio.getValue();
    m_fTPLMaxFalsePositiveRatio = oArgSwitchTPLMaxFalsePositiveRatio.getValue();
    m_fTPLMaxFalsePositiveRatioModerator = oArgSwitchTPLMaxFalsePositiveRatioModerator.getValue();

    return true;
}

/**
   @brief    This method check if command line arguments are coherent between us.
   @details  If an incoherence is detected into command line arguments a message will be displayed to user and in better cases consistency will be restored.
*/
bool CSegPerfApp::checkParamsCoherence()
{
    bool bFault = false;

    if(!(m_bSegmentationEvaluation || m_bAdvancedEvaluation || m_bSurfaceEvaluation || m_bLesionsDetectionEvaluation))
    {
        m_bSegmentationEvaluation = true;
    }

    if (m_bAdvancedEvaluation && !(m_bSegmentationEvaluation || m_bSurfaceEvaluation))
    {
        std::cout << "Switch \"advancedEvaluation\" need \"SegmentationEvaluationMetrics\" or/and SurfaceEvaluation switches." << std::endl;
        std::cout << "!!!\"SegmentationEvaluationMetrics\" has been enabled by default!!!" << std::endl;
        m_bSegmentationEvaluation = true;
    }

    if (m_fDetectionLesionMinVolume<0)
    {
        std::cout << "!!!!! Error on DetectionLesionMinVolume!!!!!" << std::endl;
        bFault = true;
    }

    if (m_fTPLMinOverlapRatio<=0 || m_fTPLMinOverlapRatio>1)
    {
        std::cout << "!!!!! Error on TPLMinOverlapRatio!!!!!" << std::endl;
        bFault = true;
    }

    if (m_fTPLMaxFalsePositiveRatio<=0 || m_fTPLMaxFalsePositiveRatio>1)
    {
        std::cout << "!!!!! Error on TPLMaxFalsePositiveRatio!!!!!" << std::endl;
        bFault = true;
    }

    if (m_fTPLMaxFalsePositiveRatioModerator<=0 || m_fTPLMaxFalsePositiveRatioModerator>1)
    {
        std::cout << "!!!!! Error on TPLMaxFalsePositiveRatioModerator!!!!!" << std::endl;
        bFault = true;
    }

    if (bFault)
    {
        std::cout << "*** 0 < DetectionLesionMinVolume               ***" << std::endl;
        std::cout << "*** 0 < TPLMinOverlapRatio                <= 1 ***" << std::endl;
        std::cout << "*** 0 < TPLMaxFalsePositiveRatio          <= 1 ***" << std::endl;
        std::cout << "*** 0 < TPLMaxFalsePositiveRatioModerator <= 1 ***" << std::endl;
    }

    return !bFault;
}

/**
   @brief    This method define an output way if none are defined by command line.
   @details  For the moment the default output way is text file.
*/
void CSegPerfApp::checkOutputCoherence()
{
    if (!(m_bTxt || m_bXml || m_bScreen))
    {
        m_bTxt = true;
    }
}

/**
   @brief    This method define the base file name for output way.
*/
void CSegPerfApp::prepareOutput()
{
    if (m_oStrBaseOut.empty())
    {
        int iDotPos = m_oStrInImage.rfind('.');
        int iFolderPos = -1;
        int iSlashPos = m_oStrInImage.rfind('/');
        int iBkSlashPos = m_oStrInImage.rfind('\\');

        if (iSlashPos != std::string::npos && iBkSlashPos != std::string::npos)
        {
            iFolderPos = iSlashPos>iBkSlashPos ? iSlashPos : iBkSlashPos;
        }
        else if (iSlashPos != std::string::npos)
        {
            iFolderPos = iSlashPos;
        }
        else if (iBkSlashPos != std::string::npos)
        {
            iFolderPos = iBkSlashPos;
        }

        ++iFolderPos;

        m_pchOutBase = "";
        m_pchOutBase.append(m_oStrInImage.begin() + iFolderPos, m_oStrInImage.end());

        if(iDotPos != std::string::npos)
            m_pchOutBase[iDotPos-iFolderPos] = 0;
    }
    else
        m_pchOutBase = m_oStrBaseOut;
}

/**
   @brief    This method execute images filter to obtain desired measures, marks and scores.
*/
void CSegPerfApp::play()
{
    long lRes = 0;

    CAnalyzer oAnalyzer(m_oStrInImage, m_oStrRefImage, m_bAdvancedEvaluation);

    if(!oAnalyzer.checkImagesMatrixAndVolumes())
        throw itk::ExceptionObject(__FILE__, __LINE__, "Orientation matrices and volumes do not match");

    if (m_iNbThreads>0)
        oAnalyzer.setNbThreads(m_iNbThreads);

    int nbLabels = oAnalyzer.getNumberOfClusters();

    std::string sOutBase = m_pchOutBase;
    std::string sOut = sOutBase;

    int i = 0;

    //////////////////////////////////////////////////////////////////////////
    // First step of loop is for global, if next steps exist it's for each layer
    do
    {
        sOut = sOutBase;
        std::stringstream outBaseTemp;

        if(i == 0)
            outBaseTemp << sOut << "_global";
        else
            outBaseTemp << sOut << "_cluster" << i;

        sOut = outBaseTemp.str();

        CResults oRes(sOut);

        processAnalyze(oAnalyzer, i);
        storeMetricsAndMarks(oRes);
        lRes = writeStoredMetricsAndMarks(oRes);

        i++;
    } while (i < nbLabels && m_bAdvancedEvaluation);
}

/**
   @brief    This method provides computing metrics, marks and scores for desired label.
   @param    [in] pi_roAnalyzer is the reference on the image analyzer class instance.
   @param    [in] pi_iIndex is label index.
   @details  This method is called by play method.
*/
void CSegPerfApp::processAnalyze(CAnalyzer&pi_roAnalyzer, int pi_iIndex)
{
    pi_roAnalyzer.selectCluster(pi_iIndex);

    // Segmentation evaluation
    if(m_bSegmentationEvaluation || m_bSurfaceEvaluation)
    {
        pi_roAnalyzer.computeITKMeasures();
        m_fSensitivity = pi_roAnalyzer.getSensitivity();
        m_fSpecificity = pi_roAnalyzer.getSpecificity();
        m_fPPV = pi_roAnalyzer.getPPV();
        m_fNPV = pi_roAnalyzer.getNPV();
        m_fDice = pi_roAnalyzer.getDiceCoefficient();
        m_fJaccard = pi_roAnalyzer.getJaccardCoefficient();
        m_fRVE = pi_roAnalyzer.getRelativeVolumeError();

        //Surfaces distances computing (Haussdorf, meanDist, Average ...)
        if(m_bSurfaceEvaluation)
        {
            m_fHausdorffDist = pi_roAnalyzer.computeHausdorffDist();
            m_fMeanDist = pi_roAnalyzer.computeMeanDist();
            m_fAverageDist = pi_roAnalyzer.computeAverageSurfaceDistance();
        }
    }

    // Detection lesions
    if(m_bLesionsDetectionEvaluation)
    {
        pi_roAnalyzer.setDetectionThresholdAlpha(m_fTPLMinOverlapRatio);
        pi_roAnalyzer.setDetectionThresholdBeta(m_fTPLMaxFalsePositiveRatio);
        pi_roAnalyzer.setDetectionThresholdGamma(m_fTPLMaxFalsePositiveRatioModerator);
        pi_roAnalyzer.setMinLesionVolumeDetection(m_fDetectionLesionMinVolume);
        pi_roAnalyzer.getDetectionMarks(m_fPPVL, m_fSensL, m_fF1);
    }
}

/**
   @brief    This method store results into CResults class instance.
   @param    [in] pi_roRes is the reference output writer class instance.
   @details  This method is called by play method after processAnalyze method.
*/
void CSegPerfApp::storeMetricsAndMarks(CResults&pi_roRes)
{
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // Following code put out results
    //////////////////////////////////////////////////////////////////////////
    pi_roRes.setTxt(m_bTxt);
    pi_roRes.setXml(m_bXml);
    pi_roRes.setScreen(m_bScreen);

    //Segmentation
    if(m_bSegmentationEvaluation)
    {
        pi_roRes.activeMeasurementOutput(CResults::eMesureDice);
        pi_roRes.setDice(m_fDice);
        pi_roRes.activeMeasurementOutput(CResults::eMesureJaccard);
        pi_roRes.setJaccard(m_fJaccard);
        pi_roRes.activeMeasurementOutput(CResults::eMesureSensibility);
        pi_roRes.setSensibility(m_fSensitivity);
        pi_roRes.activeMeasurementOutput(CResults::eMesureSpecificity);
        pi_roRes.setSpecificity(m_fSpecificity);
        pi_roRes.activeMeasurementOutput(CResults::eMesureNPV);
        pi_roRes.setNPV(m_fNPV);
        pi_roRes.activeMeasurementOutput(CResults::eMesurePPV);
        pi_roRes.setPPV(m_fPPV);
        pi_roRes.activeMeasurementOutput(CResults::eMesureRelativeVolumeError);
        pi_roRes.setRVE(m_fRVE*100);
    }

    //Surfaces distances
    if(m_bSurfaceEvaluation)
    {
        pi_roRes.activeMeasurementOutput(CResults::eMesureDistHausdorff);
        pi_roRes.setHausdorffDist(m_fHausdorffDist);
        pi_roRes.activeMeasurementOutput(CResults::eMesureDistMean);
        pi_roRes.setContourMeanDist(m_fMeanDist);
        pi_roRes.activeMeasurementOutput(CResults::eMesureDistAverage);
        pi_roRes.setAverageSurfaceDist(m_fAverageDist);
    }

    //Lesion detection
    if(m_bLesionsDetectionEvaluation)
    {
        pi_roRes.activeMeasurementOutput(CResults::eMesurePPVL);
        pi_roRes.setPPVL(m_fPPVL);
        pi_roRes.activeMeasurementOutput(CResults::eMesureSensL);
        pi_roRes.setSensL(m_fSensL);
        pi_roRes.activeMeasurementOutput(CResults::eMesureF1Test);
        pi_roRes.setF1test(m_fF1);
    }
}

/**
   @brief    This method flush CResults class instance.
   @param    [in] pi_roRes is the reference output writer class instance.
   @details  This method is called by play method after storeMetricsAndMarks method.
*/
long CSegPerfApp::writeStoredMetricsAndMarks(CResults&pi_roRes)
{
    long lRes = (long) !pi_roRes.save();
    return lRes;
}

/**
   @brief    This method display information about SegPerfAnalyzer results.
*/
void CSegPerfApp::about()
{
    std::cout << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "SegPerfAnalyser (Segmentation Performance Analyzer) provides different" << std::endl;
    std::cout << "marks, metrics and scores for segmentation evaluation." << std::endl;
    std::cout << std::endl;
    std::cout << "3 categories are available:" << std::endl;
    std::cout << "    - SEGMENTATION EVALUATION:" << std::endl;
    std::cout << "        Dice, the mean overlap" << std::endl;
    std::cout << "        Jaccard, the union overlap" << std::endl;
    std::cout << "        Sensitivity" << std::endl;
    std::cout << "        Specificity" << std::endl;
    std::cout << "        NPV (Negative Predictive Value)" << std::endl;
    std::cout << "        PPV (Positive Predictive Value)" << std::endl;
    std::cout << "        RVE (Relative Volume Error) in percentage" << std::endl;
    std::cout << "    - SURFACE DISTANCE EVALUATION:" << std::endl;
    std::cout << "        Hausdorff distance" << std::endl;
    std::cout << "        Contour mean distance" << std::endl;
    std::cout << "        Average surface distance" << std::endl;
    std::cout << "    - DETECTION LESIONS EVALUATION:" << std::endl;
    std::cout << "        PPVL (Positive Predictive Value for Lesions)" << std::endl;
    std::cout << "        SensL, Lesion detection sensitivity" << std::endl;
    std::cout << "        F1 Score, a F1 Score between PPVL and SensL" << std::endl;
    std::cout << std::endl;


    std::cout << "Results are provided as follows: " << std::endl;
    char const*const*const ppchNameTab = CResults::getMeasureNameTable();
    for (int i=0;i < CResults::eMesureLast;++i)
        std::cout << ppchNameTab[i]<<";\t";

    std::cout << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "********************************************************************************" << std::endl;
}

