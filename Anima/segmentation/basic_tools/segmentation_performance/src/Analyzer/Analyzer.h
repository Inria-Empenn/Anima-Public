#ifndef __ANALYZER_H
#define __ANALYZER_H

#include <math.h>
#include <vector>
#include <algorithm>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include "SegmentationMeasuresImageFilter.h"
#include <itkLabelContourImageFilter.h>
#include <itkBinaryContourImageFilter.h>
#include <itkSimpleFilterWatcher.h>
#include <itkHausdorffDistanceImageFilter.h>
#include <itkDirectedHausdorffDistanceImageFilter.h>
#include <itkContourMeanDistanceImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkImageDuplicator.h>



/**
* @file Analyzer.h
* @brief Description of CAnalyzer class. The class to compute various metrics to evaluate segmentation results.
* @author Baptiste LAURENT
* @author Florent  LERAY
* @date 12/04/2016
* @version 2.0
*/



/**
* @class CAnalyzer Analyzer.h "Class to compute various metrics to evaluate segmentation results."
* @brief Class to compute various metrics to evaluate segmentation results.
*/
class CAnalyzer
{
   //////////////////////////////////////////////////////////////////////////
   //  Les membres
private:
   unsigned int m_uiNbLabels;   /*!<Number of Labels. */
   bool m_bValuesComputed;      /*!<Boolean to check if values have been computed. */
   bool m_bContourDetected;     /*!<Boolean to check if contour detection have been done. */


   double m_dfDetectionThresholdAlpha;
   double m_dfDetectionThresholdBeta;
   double m_dfDetectionThresholdGamma;
   double m_dfMinLesionVolumeDetection;


   typedef itk::Image <unsigned short, 3> ImageType;
   typedef itk::ImageFileReader <ImageType> ImageReaderType;
   typedef itk::ImageRegionConstIterator <ImageType> ImageIteratorType;
   typedef itk::SegmentationMeasuresImageFilter<ImageType> FilterType;

   ImageType::Pointer m_imageTest;
   ImageType::Pointer m_imageRef;
   ImageType::Pointer m_imageTestContour;
   ImageType::Pointer m_imageRefContour;
   ImageType::Pointer m_imageRefDuplicated;
   ImageType::Pointer m_imageTestDuplicated;

   FilterType::Pointer m_oFilter;
   itk::ThreadIdType m_ThreadNb;





protected:

public:


   //////////////////////////////////////////////////////////////////////////
   //  Les methodes
public:
   CAnalyzer(char *pi_pchImageTestName, char *pi_pchImageRefName, bool advancedEvaluation);
   ~CAnalyzer();

   bool checkImagesMatixAndVolumes();

   void setNbThreads(int pi_iNbThreads);
   void selectCluster(unsigned int);

   void setDetectionThresholdAlpha(double pi_fVal)
   {
      m_dfDetectionThresholdAlpha = pi_fVal;
   }

   void setDetectionThresholdBeta(double pi_fVal)
   {
      m_dfDetectionThresholdBeta = pi_fVal;
   }

   void setDetectionThresholdGamma(double pi_fVal)
   {
      m_dfDetectionThresholdGamma = pi_fVal;
   }

   void setMinLesionVolumeDetection(double pi_fVal)
   {
      m_dfMinLesionVolumeDetection = pi_fVal;
   }

   float computeHausdorffDist();
   float computeMeanDist();
   float computeAverageSurfaceDistance();
   void  computeITKMeasures();

   //getters
   float getUnionOverlap();
   float getMeanOverlap();
   float getSensitivity();
   float getSpecificity();
   float getPPV();
   float getNPV();
   float getDiceCoefficient();
   float getJaccardCoefficient();
   float getRelativeVolumeError();
   bool getDetectionMarks(float&po_fPPVL, float&po_fSensL, float&po_fF1);



   int getNumberOfClusters();

protected:
   void formatLabels();
   void contourDectection();
   void checkNumberOfLabels(int, int);

   int getTruePositiveLesions(int pi_iNbLabelsRef, int pi_iNbLabelsTest, int * *pi_ppiOverlapTab);
   bool falsePositiveRatioTester(int pi_iLessionReference, int pi_iNbLabelsTest, int pi_iNbLabelsRef, int * *pi_ppiOverlapTab, int *pi_piTPRowSumTab, int *pi_piColumnSumTab, double pi_dfBeta, double pi_dfGamma);
   void getOverlapTab(int&po_iNbLabelsRef, int&po_iNbLabelsTest, int* *&po_ppiOverlapTab);

private:
   CAnalyzer();
   void transposer(int pi_iNbLabelsRef, int pi_iNbLabelsTest, int * *pi_ppiOverlapTab, int* *&po_rppiOverlapTabTransposed);
   void removeOverlapTab(int * *pi_ppiOverlapTab, int pi_iNbLabelsRef);
};



#endif