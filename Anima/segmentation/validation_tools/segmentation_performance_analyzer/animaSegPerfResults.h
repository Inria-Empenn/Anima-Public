#pragma once
#include <string>

namespace anima
{

/**
* @class SegPerfResults
* @brief Class to format and saves results.
*/
class SegPerfResults
{
public:
    typedef enum
    {
        eMesureJaccard = 0,
        eMesureDice,
        eMesureSensibility,
        eMesureSpecificity,
        eMesurePPV,
        eMesureNPV,
        eMesureRelativeVolumeError,
        eMesureDistHausdorff,
        eMesureDistMean,
        eMesureDistAverage,
        eMesurePPVL,
        eMesureSensL,
        eMesureF1Test,
        eMesureLast
    }eMesureName;

    SegPerfResults(std::string &pi_pchBaseFileName);
    ~SegPerfResults();

    bool save();

    bool activeMeasurementOutput(eMesureName pi_eVal);

    static char const*const*const getMeasureNameTable();

    /**
      @brief    Enable or disable text file results.
      @param    [in] pi_bEnable Enable or disable.
   */
    void setTxt(bool pi_bEnable = true)
    {
        m_bTxt = pi_bEnable;
    }

    /**
      @brief    Enable or disable XML file results.
      @param    [in] pi_bEnable Enable or disable.
   */
    void setXml(bool pi_bEnable = true)
    {
        m_bXml = pi_bEnable;
    }

    /**
   @brief    Enable or disable on screen results.
   @param    [in] pi_bEnable Enable or disable.
   */
    void setScreen(bool pi_bEnable = true)
    {
        m_bScreen = pi_bEnable;
    }

    /**
      @brief    Set the result value of Jaccard measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setJaccard(double pi_fVal)
    {
        m_fResTab[eMesureJaccard] = pi_fVal;
    }

    /**
      @brief    Set the result value of Dice measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setDice(double pi_fVal)
    {
        m_fResTab[eMesureDice] = pi_fVal;
    }

    /**
      @brief    Set the result value of Sensibility measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setSensibility(double pi_fVal)
    {
        m_fResTab[eMesureSensibility] = pi_fVal;
    }

    /**
      @brief    Set the result value of Specificity measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setSpecificity(double pi_fVal)
    {
        m_fResTab[eMesureSpecificity] = pi_fVal;
    }

    /**
      @brief    Set the result value of PPV (Positive Predictive Value) measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setPPV(double pi_fVal)
    {
        m_fResTab[eMesurePPV] = pi_fVal;
    }

    /**
      @brief    Set the result value of NPV (Negative Predictive Value) measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setNPV(double pi_fVal)
    {
        m_fResTab[eMesureNPV] = pi_fVal;
    }

    /**
      @brief    Set the result value of Relative volume error.
      @param    [in] pi_fVal Measure result value.
   */
    void setRVE(double pi_fVal)
    {
        m_fResTab[eMesureRelativeVolumeError] = pi_fVal;
    }

    /**
      @brief    Set the result value of DistHausdorff measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setHausdorffDist(double pi_fVal)
    {
        m_fResTab[eMesureDistHausdorff] = pi_fVal;
    }

    /**
      @brief    Set the result value of contour mean distance measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setContourMeanDist(double pi_fVal)
    {
        m_fResTab[eMesureDistMean] = pi_fVal;
    }

    /**
      @brief    Set the result value of average surface distance measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setAverageSurfaceDist(double pi_fVal)
    {
        m_fResTab[eMesureDistAverage] = pi_fVal;
    }

    /**
   @brief    Set the result value of PPVL measure.
   @param    [in] pi_fVal Measure result value.
*/
    void setPPVL(double pi_fVal)
    {
        m_fResTab[eMesurePPVL] = pi_fVal;
    }

    /**
      @brief    Set the result value of SensL measure.
      @param    [in] pi_fVal Measure result value.
   */
    void setSensL(double pi_fVal)
    {
        m_fResTab[eMesureSensL] = pi_fVal;
    }

    /**
      @brief    Set the result value of F1 score of F-test.
      @param    [in] pi_fVal Measure result value.
   */
    void setF1test(double pi_fVal)
    {
        m_fResTab[eMesureF1Test] = pi_fVal;
    }

private:
    SegPerfResults();

    double m_fResTab[eMesureLast];       /*!<Table of measurements results. */
    bool m_bResActiveTab[eMesureLast];  /*!<Table of measurements present in outputs*/
    bool m_bTxt;                        /*!<Enable txt output results format. */
    bool m_bXml;                        /*!<Enable Xml output results format. */
    bool m_bScreen;                     /*!<Enable screen output results format. */
    std::string m_pchBaseOutputFileName; /*!<Base name for output results file. */

    static char const*const m_ppchMeasureNameTable[eMesureLast];   /*!<Table to associate a measure name for each measure index defined into eMesureName */
};

} // end namespace anima
