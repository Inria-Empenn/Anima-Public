#include "Analyzer.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <exception>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageIterator.h>
#include <itkMultiThreader.h>
#include <itkImageDuplicator.h>

/**
   @brief    Constructor.
   @details  Read input images
   @param	[in] pi_pchImageTestName Name of the image to evaluate
   @param	[in] pi_pchImageRefName Name of the ground truth image
   @param	[in] bAdvancedEvaluation
*/
CAnalyzer::CAnalyzer(std::string &pi_pchImageTestName, std::string &pi_pchImageRefName, bool bAdvancedEvaluation)
{
    m_ThreadNb = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();

    m_dfDetectionThresholdAlpha = 0.05;
    m_dfDetectionThresholdBeta = 0.50;
    m_dfDetectionThresholdGamma = 0.50;
    m_dfMinLesionVolumeDetection = 3.0;

    m_uiNbLabels = 0;

    ImageReaderType::Pointer imageTestReader = ImageReaderType::New();
    imageTestReader->SetFileName(pi_pchImageTestName);

    try
    {
        imageTestReader->Update();
    }

    catch (itk::ExceptionObject& e)
    {
        std::cerr << "exception in file reader " << std::endl;
        std::cerr << e << std::endl;
        return;
    }

    m_imageTest = imageTestReader->GetOutput();

    ImageReaderType::Pointer imageRefReader = ImageReaderType::New();
    imageRefReader->SetFileName(pi_pchImageRefName);

    try
    {
        imageRefReader->Update();
    }

    catch (itk::ExceptionObject& e)
    {
        std::cerr << "exception in file reader " << std::endl;
        std::cerr << e << std::endl;
        return;
    }

    m_imageRef = imageRefReader->GetOutput();

    if (!checkImagesMatrixAndVolumes())
        throw std::runtime_error("Images are incompatible");

    this->formatLabels();

    if(bAdvancedEvaluation && m_uiNbLabels > 2)
    {
        //TO DO : attention si pas meme nombre de labels
        typedef itk::ImageDuplicator< ImageType > DuplicatorType;
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(m_imageRef);
        duplicator->Update();
        m_imageRefDuplicated = ImageType::New();
        m_imageRefDuplicated = duplicator->GetOutput();

        duplicator->SetInputImage(m_imageTest);
        duplicator->Update();
        m_imageTestDuplicated = ImageType::New();
        m_imageTestDuplicated = duplicator->GetOutput();
    }

    m_bValuesComputed = false;
    m_bContourDetected = false;
}

/**
@brief Check if the 2 inputs images are compatible.
@return	true if image are compatible else false.
*/
bool CAnalyzer::checkImagesMatrixAndVolumes()
{
    bool bRes = true;

    unsigned int uiTestedImageDimension = m_imageTest->GetImageDimension();
    const itk::Vector<ImageType::SpacingValueType, ImageType::ImageDimension> oVectSpacingTested = m_imageTest->GetSpacing();
    const itk::Size<ImageType::ImageDimension> oSizeTested = m_imageTest->GetLargestPossibleRegion().GetSize();
    const itk::Point<ImageType::IndexValueType, ImageType::ImageDimension> oPointOriginTested = m_imageTest->GetOrigin();

    unsigned int uiRefImageDimension = m_imageRef->GetImageDimension();
    const itk::Vector<ImageType::SpacingValueType, ImageType::ImageDimension> oVectSpacingRef = m_imageRef->GetSpacing();
    const itk::Size<ImageType::ImageDimension> oSizeRef = m_imageRef->GetLargestPossibleRegion().GetSize();
    const itk::Point<ImageType::IndexValueType, ImageType::ImageDimension> oPointOriginRef = m_imageRef->GetOrigin();

    bRes = uiTestedImageDimension == uiRefImageDimension;

    if (bRes)
    {
        //Check dimensions
        for (int i = 0; i < uiTestedImageDimension; ++i)
        {
            if (oSizeTested[i] != oSizeRef[i])
            {
                bRes = false;
                std::cerr << "Different image dimensions " << i << " Tested size : " << oSizeTested[i] << " Reference size : " << oSizeRef[i] << std::endl;
            }
        }

        //Check origin
        for (int i = 0; i < uiTestedImageDimension; ++i)
        {
            if (oPointOriginTested[i] != oPointOriginRef[i])
            {
                bRes = false;
                std::cerr << "Different image origin " << i << " Tested origin : " << oPointOriginTested[i] << " Reference origin : " << oPointOriginRef[i] << std::endl;
            }
        }

        //Check direction
        for (int i = 0; i < uiTestedImageDimension; ++i)
        {
            for (int j = 0; j < uiTestedImageDimension; ++j)
            {
                double dDiff = m_imageTest->GetDirection().GetVnlMatrix().get(i, j) - m_imageRef->GetDirection().GetVnlMatrix().get(i, j);


                if ((dDiff > 0.000001) || (dDiff < -0.000001))
                {
                    bRes = false;
                    std::cerr << "Different image direction " << i << " - " << j << std::endl;
                    std::cerr << "Diff = " << dDiff << " Tested direction : " << m_imageTest->GetDirection().GetVnlMatrix().get(i, j) << " Reference direction : " << m_imageRef->GetDirection().GetVnlMatrix().get(i, j) << std::endl;
                }
            }
        }
    }
    else
    {
        std::cerr << "Different image dimensions." << std::endl << "Tested is in  : " << uiTestedImageDimension << "D and Ref is in :" << uiRefImageDimension << "D" << std::endl;
    }

    return bRes;
}

CAnalyzer::~CAnalyzer()
{
}

/**
   @brief    Return the number of clusters.
   @return	Number of clusters
*/
int CAnalyzer::getNumberOfClusters()
{
    return this->m_uiNbLabels;
}

/**
   @brief    Select the cluster we want to use to compute evaluation results.
   @param	[in] iCluster
*/
void CAnalyzer::selectCluster(unsigned int iCluster)
{
    //cout << "select cluster : "<< iCluster<<endl;
    if(iCluster!=0)
    {
        typedef itk::ImageDuplicator< ImageType > DuplicatorType;
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(m_imageRefDuplicated);
        duplicator->Update();
        m_imageRef = duplicator->GetOutput();

        duplicator->SetInputImage(m_imageTestDuplicated);
        duplicator->Update();
        m_imageTest = duplicator->GetOutput();

        unsigned int dimX = m_imageRef->GetLargestPossibleRegion().GetSize()[0];
        unsigned int dimY = m_imageRef->GetLargestPossibleRegion().GetSize()[1];
        unsigned int dimZ = m_imageRef->GetLargestPossibleRegion().GetSize()[2];

        ImageType::IndexType index;

        for(int z=0; z<dimZ; z++)
        {
            for(int y=0; y<dimY; y++)
            {
                for(int x=0; x<dimX; x++)
                {
                    index[0]=x;
                    index[1]=y;
                    index[2]=z;

                    int valueRef = m_imageRef->GetPixel(index);
                    int valueTest = m_imageTest->GetPixel(index);

                    if(valueRef != iCluster)
                        m_imageRef->SetPixel(index, 0);
                    else
                        m_imageRef->SetPixel(index, 1);

                    if(valueTest != iCluster)
                        m_imageTest->SetPixel(index, 0);
                    else
                        m_imageTest->SetPixel(index, 1);
                }
            }
        }
        this->m_uiNbLabels = 2;
    }

    return;
}

/**
   @brief    Extract the contour of the image to evaluate and the ground truth
*/
void CAnalyzer::contourDectection()
{
    //cout<<"****** CONTOUR DETECTION ******"<<endl;
    //cout<<"nbLabels = "<<this->m_uiNbLabels<<endl;
    if(this->m_uiNbLabels == 2)
    {
        typedef itk::BinaryContourImageFilter< ImageType, ImageType > FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput( m_imageTest );
        filter->SetFullyConnected( 0 );
        filter->SetForegroundValue( 1 );
        filter->SetBackgroundValue( 0 );
        filter->Update();


        m_imageTestContour = ImageType::New();
        m_imageTestContour = filter->GetOutput();

        FilterType::Pointer filter2 = FilterType::New();
        filter2->SetInput( m_imageRef );
        filter2->SetFullyConnected( 0 );
        filter2->SetForegroundValue( 1 );
        filter2->SetBackgroundValue( 0 );
        filter2->Update();

        m_imageRefContour = ImageType::New();
        m_imageRefContour = filter2->GetOutput();
    }

    else
    {
        typedef itk::LabelContourImageFilter< ImageType, ImageType > FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput( m_imageTest );
        filter->SetFullyConnected( 0 );
        filter->SetBackgroundValue( 0 );
        filter->Update();

        m_imageTestContour = ImageType::New();
        m_imageTestContour = filter->GetOutput();

        FilterType::Pointer filter2 = FilterType::New();
        filter2->SetInput( m_imageRef );
        filter2->SetFullyConnected( 0 );
        filter2->SetBackgroundValue( 0 );
        filter2->Update();

        m_imageRefContour = ImageType::New();
        m_imageRefContour = filter2->GetOutput();
    }
    m_bContourDetected = true;
}

/**
   @brief    Compute different measures with ITK to evaluate segmentation
*/
void CAnalyzer::computeITKMeasures()
{
#ifdef _DEBUG
    cout<<endl<<"****** COMPUTE ITK MEASURES ******"<<endl;
#endif

    m_oFilter = FilterType::New();

    m_oFilter->SetNumberOfThreads(m_ThreadNb);
    m_oFilter->SetSourceImage( m_imageTest);
    m_oFilter->SetTargetImage( m_imageRef );

#ifdef _DEBUG
    std::cout<<"start compute itk measures"<<std::endl;
#endif

    try
    {
        m_oFilter->Update();
    }
    catch( itk::ExceptionObject& excp  )
    {
        std::cerr << excp << std::endl;
        return;
    }

#ifdef _DEBUG
    std::cout<<"end compute itk measures"<<std::endl;
    std::cout<< "UnionOverlap (Jaccard) : " <<m_oFilter->getUnionOverlap()<<std::endl;
    std::cout<< "MeanOverlap (Dice)     : " <<m_oFilter->getMeanOverlap()<<std::endl;
    cout<<"****** COMPUTE ITK MEASURES end ******"<<endl;
#endif
}

/**
   @brief  Getter of Union overlap
   @return Union overlap in float
*/
float CAnalyzer::getUnionOverlap()
{
    return m_oFilter->getUnionOverlap();
}

/**
   @brief  Getter of Mean overlap
   @return Mean overlap in float
*/
float CAnalyzer::getMeanOverlap()
{
    return m_oFilter->getMeanOverlap();
}

/**
   @brief  Getter of Sensibility
   @return Sensibility in float
*/
float CAnalyzer::getSensitivity()
{
    return m_oFilter->getSensitivity();
}

/**
   @brief  Getter of Specificity
   @return Specificity in float
*/
float CAnalyzer::getSpecificity()
{
    return m_oFilter->getSpecificity();
}

/**
   @brief  Getter of Positive predictive value
   @return Positive predictive value in float
*/
float CAnalyzer::getPPV()
{
    return m_oFilter->getPPV();
}

/**
   @brief  Getter of Negative predictive value
   @return Negative predictive value in float
*/
float CAnalyzer::getNPV()
{
    return m_oFilter->getNPV();
}

/**
   @brief  Getter of Dice coefficient
   @return Dice coefficient in float
*/
float CAnalyzer::getDiceCoefficient()
{
    return m_oFilter->GetDiceCoefficient();
}

/**
   @brief  Getter of Jaccard coefficient
   @return Jaccard coefficient in float
*/
float CAnalyzer::getJaccardCoefficient()
{
    return m_oFilter->GetJaccardCoefficient();
}

/**
   @brief  Getter of Relative volume error
   @return Relative volume error in float
*/
float CAnalyzer::getRelativeVolumeError()
{
    return m_oFilter->getRelativeVolumeError();
}

/**
   @brief    Check if the number of labels is the same for both input images
*/
void CAnalyzer::checkNumberOfLabels(int iNbLabelsImageTest, int iNbLabelsImageRef)
{
    if(iNbLabelsImageTest<=1)
    {
        m_uiNbLabels = 1;
        std::cerr << "ERROR : Number of labels for ground truth is 0 !" << std::endl;
        return;
    }

    if(iNbLabelsImageRef<=1)
    {
        m_uiNbLabels = 1;
        std::cerr << "ERROR : Number of labels for reference image is 0 !"<< std::endl;
        return;
    }

    if (iNbLabelsImageTest == iNbLabelsImageRef)
    {
        m_uiNbLabels = iNbLabelsImageTest;
    }
    else
    {
        if(iNbLabelsImageTest > 2)
        {
            //binariser ImageTest
            unsigned int dimX = m_imageTest->GetLargestPossibleRegion().GetSize()[0];
            unsigned int dimY = m_imageTest->GetLargestPossibleRegion().GetSize()[1];
            unsigned int dimZ = m_imageTest->GetLargestPossibleRegion().GetSize()[2];

            ImageType::IndexType index;

            for(int z=0; z<dimZ; z++)
            {
                for(int y=0; y<dimY; y++)
                {
                    for(int x=0; x<dimX; x++)
                    {
                        index[0]=x;
                        index[1]=y;
                        index[2]=z;

                        if( m_imageTest->GetPixel(index) != 0 )
                            m_imageTest->SetPixel(index, 1);
                    }
                }
            }

            m_uiNbLabels = 2;
            std::cout << "WARNING : Segmented image have not the same number of labels as ground truth, it have been binarized for segmentation evaluation" << std::endl;
        }

        if(iNbLabelsImageRef > 2)
        {
            //binariser ImageRef

            unsigned int dimX = m_imageRef->GetLargestPossibleRegion().GetSize()[0];
            unsigned int dimY = m_imageRef->GetLargestPossibleRegion().GetSize()[1];
            unsigned int dimZ = m_imageRef->GetLargestPossibleRegion().GetSize()[2];

            ImageType::IndexType index;

            for(int z=0; z<dimZ; z++)
            {
                for(int y=0; y<dimY; y++)
                {
                    for(int x=0; x<dimX; x++)
                    {
                        index[0]=x;
                        index[1]=y;
                        index[2]=z;

                        if( m_imageRef->GetPixel(index) != 0 )
                            m_imageRef->SetPixel(index, 1);
                    }
                }
            }

            m_uiNbLabels = 2;
            std::cout << "WARNING : Ground truth have not the same number of labels as segmented image, both images have been binarized for segmentation evaluation" << std::endl;
        }
    }
}

/**
   @brief    Format labels values for the image to evaluate and the ground truth to obtain : background pixels values = 0, cluster 1 pixels values = 1 ...etc
*/
void CAnalyzer::formatLabels()
{
    ImageIteratorType refIt (m_imageTest, m_imageTest->GetLargestPossibleRegion());

    std::vector <unsigned int> usefulLabels;

    while (!refIt.IsAtEnd())
    {
        bool isAlreadyIn = false;
        for (unsigned int i = 0; i < usefulLabels.size(); ++i)
        {
            if (refIt.Get() == usefulLabels[i])
            {
                isAlreadyIn = true;
                break;
            }
        }

        if (!isAlreadyIn)
        {
            usefulLabels.push_back(refIt.Get());
        }

        ++refIt;
    }

    std::sort(usefulLabels.begin(), usefulLabels.end());

    unsigned int dimX = m_imageTest->GetLargestPossibleRegion().GetSize()[0];
    unsigned int dimY = m_imageTest->GetLargestPossibleRegion().GetSize()[1];
    unsigned int dimZ = m_imageTest->GetLargestPossibleRegion().GetSize()[2];

    ImageType::IndexType index;

    for(int z=0; z<dimZ; z++)
    {
        for(int y=0; y<dimY; y++)
        {
            for(int x=0; x<dimX; x++)
            {
                index[0]=x;
                index[1]=y;
                index[2]=z;

                int indexvalue = std::find(usefulLabels.begin(), usefulLabels.end(), m_imageTest->GetPixel(index))-usefulLabels.begin();
                m_imageTest->SetPixel(index, indexvalue);
            }
        }
    }


    //imageRef

    ImageIteratorType refItRef(m_imageRef, m_imageRef->GetLargestPossibleRegion());

    std::vector <unsigned int> usefulLabels2;

    while (!refItRef.IsAtEnd())
    {
        bool isAlreadyIn = false;
        for (unsigned int i = 0; i < usefulLabels2.size(); ++i)
        {
            if (refItRef.Get() == usefulLabels2[i])
            {
                isAlreadyIn = true;
                break;
            }
        }

        if (!isAlreadyIn) {
            usefulLabels2.push_back(refItRef.Get());
        }

        ++refItRef;
    }

    std::sort(usefulLabels2.begin(), usefulLabels2.end());



    dimX = m_imageRef->GetLargestPossibleRegion().GetSize()[0];
    dimY = m_imageRef->GetLargestPossibleRegion().GetSize()[1];
    dimZ = m_imageRef->GetLargestPossibleRegion().GetSize()[2];


    for(int z=0; z<dimZ; z++)
    {
        for(int y=0; y<dimY; y++)
        {
            for(int x=0; x<dimX; x++)
            {
                index[0]=x;
                index[1]=y;
                index[2]=z;

                int indexvalue = std::find(usefulLabels2.begin(), usefulLabels2.end(), m_imageRef->GetPixel(index))-usefulLabels2.begin();
                m_imageRef->SetPixel(index, indexvalue);
            }
        }
    }

    int iNbOfLabelsImageTest = usefulLabels.size();
    int iNbOfLabelsImageRef = usefulLabels2.size();

    checkNumberOfLabels(iNbOfLabelsImageTest, iNbOfLabelsImageRef);

    return;
}

/**
@brief  Compute Haussdorf distance
@return hausdorffDistance in float
*/
float CAnalyzer::computeHausdorffDist()
{
    FilterType::RealType hausdorffDistance = std::numeric_limits<float>::quiet_NaN();

    if (m_uiNbLabels>1)
    {
        // compute the Hausdorff distance
        typedef itk::HausdorffDistanceImageFilter<ImageType, ImageType> FilterType;
        FilterType::Pointer filter = FilterType::New();

        filter->SetInput1(m_imageTest);
        filter->SetInput2(m_imageRef);

        try
        {
            filter->Update();
        }

        catch (itk::ExceptionObject& e)
        {
            std::cerr << "exception in filter " << std::endl;
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
        }

        hausdorffDistance = filter->GetHausdorffDistance();
    }

    return hausdorffDistance;
}

/**
@brief   Compute mean distance
@return  meanDistance
*/
float CAnalyzer::computeMeanDist()
{
    FilterType::RealType meanDistance = std::numeric_limits<float>::quiet_NaN();

    if (m_uiNbLabels > 1)
    {
        if (!this->m_bContourDetected)
            this->contourDectection();

        // compute the Hausdorff distance H(image1,image2)
        typedef itk::ContourMeanDistanceImageFilter<ImageType, ImageType> FilterType;
        FilterType::Pointer filter = FilterType::New();

        filter->SetInput1(m_imageTestContour);
        filter->SetInput2(m_imageRefContour);

        try
        {
            filter->Update();
        }

        catch (itk::ExceptionObject& e)
        {
            std::cerr << "exception in filter " << std::endl;
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
        }

        meanDistance = filter->GetMeanDistance();
    }

    return meanDistance;
}

/**
@brief   Compute average surface distance
@return  average surface distance
*/
float CAnalyzer::computeAverageSurfaceDistance()
{
    float meanDistance = std::numeric_limits<float>::quiet_NaN();

    if (m_uiNbLabels > 1)
    {
        if (!this->m_bContourDetected)
            this->contourDectection();

        float sum_dist = 0;
        float sum_size = 0;

        for (int i = 1; i < m_uiNbLabels; i++)
        {
            std::vector <ImageType::PointType> coordRef;
            coordRef.clear();
            ImageIteratorType refContourIt(m_imageRefContour, m_imageRefContour->GetLargestPossibleRegion());

            while (!refContourIt.IsAtEnd())
            {
                if (refContourIt.Get() == i)
                {
                    ImageType::IndexType oIndex = refContourIt.GetIndex();
                    ImageType::PointType oPoint;
                    m_imageRefContour->TransformIndexToPhysicalPoint(oIndex, oPoint);
                    coordRef.push_back(oPoint);
                }
                ++refContourIt;
            }

            std::vector <ImageType::PointType> coordTest;
            coordTest.clear();

            ImageIteratorType testContourIt(m_imageTestContour, m_imageTestContour->GetLargestPossibleRegion());

            while (!testContourIt.IsAtEnd())
            {
                if (testContourIt.Get() == i)
                {
                    ImageType::IndexType oIndex = testContourIt.GetIndex();
                    ImageType::PointType oPoint;
                    m_imageTestContour->TransformIndexToPhysicalPoint(oIndex, oPoint);
                    coordTest.push_back(oPoint);
                }
                ++testContourIt;
            }

            std::vector <float> distance1, distance2;
            float fDistanceTemp1 = 1000000, fDistanceTemp2 = 1000000;
            float distanceValue = 0;
            for (int m = 0; m < coordRef.size(); m++)
            {
                for (int n = 0; n < coordTest.size(); n++)
                {
                    distanceValue = sqrt(pow((float)(coordRef[m][0] - coordTest[n][0]), 2) + pow((float)(coordRef[m][1] - coordTest[n][1]), 2) + pow((float)(coordRef[m][2] - coordTest[n][2]), 2));
                    if (distanceValue < fDistanceTemp1)
                        fDistanceTemp1 = distanceValue;
                }

                distance1.push_back(fDistanceTemp1);
            }

            for (int m = 0; m < coordTest.size(); m++)
            {
                for (int n = 0; n < coordRef.size(); n++)
                {
                    distanceValue = sqrt(pow((float)(coordTest[m][0] - coordRef[n][0]), 2) + pow((float)(coordTest[m][1] - coordRef[n][1]), 2) + pow((float)(coordTest[m][2] - coordRef[n][2]), 2));

                    if (distanceValue < fDistanceTemp2)
                        fDistanceTemp2 = distanceValue;
                }

                distance2.push_back(fDistanceTemp2);
            }

            sum_dist = sum_dist + (float)accumulate(distance1.begin(), distance1.end(), 0.0) + (float)accumulate(distance2.begin(), distance2.end(), 0.0);
            sum_size = sum_size + (float)distance1.size() + (float)distance2.size();
        }

        meanDistance = sum_dist / sum_size;
    }
    return meanDistance;
}

/**
@brief  Compute useful variables to get detection scores.
@param  [out] po_fFNL the output FNL (False Negative Lesions) ratio.
@param  [out] po_fFPL the output FPL (False Positive Lesions) ratio.
@param  [out] po_fPPVL the output PPVL (Positive Predictive Value for Lesions)  ratio.
@param  [out] po_fSensL the output SensL, Lesion detection sensitivity ratio.
@param  [out] po_fF1 the output F1_score.
@return true
*/
bool CAnalyzer::getDetectionMarks(float&po_fPPVL, float&po_fSensL, float&po_fF1)
{
    bool bRes = true;

    /*if (m_uiNbLabels>1)
   {*/
    int iNbLabelsRef = 0;
    int iNbLabelsTest = 0;
    int * *ppiOverlapTab = NULL;
    int * *ppiOverlapTabTransposed = NULL;
    int iTPLgt = 0;
    int iTPLd = 0;

    getOverlapTab(iNbLabelsRef, iNbLabelsTest, ppiOverlapTab);
    transposer(iNbLabelsRef, iNbLabelsTest, ppiOverlapTab, ppiOverlapTabTransposed);
    if (iNbLabelsRef>1 && iNbLabelsTest>1)
    {
        iTPLgt = getTruePositiveLesions(iNbLabelsRef, iNbLabelsTest, ppiOverlapTab);
        iTPLd  = getTruePositiveLesions(iNbLabelsTest, iNbLabelsRef, ppiOverlapTabTransposed);

        //Deallocates tables allocated by getOverlapTab and transposer
        removeOverlapTab(ppiOverlapTab, iNbLabelsRef);
        removeOverlapTab(ppiOverlapTabTransposed, iNbLabelsTest);
    }

    po_fPPVL = (float)((double)iTPLd / (double)(iNbLabelsTest-1));     //po_fTPLd  = (float)((double)iTPLd  / (double)(iNbLabelsTest-1));// The "-1" is to reject background label
    po_fSensL = (float)((double)iTPLgt / (double)(iNbLabelsRef-1));    //po_fTPLgt = (float)((double)iTPLgt / (double)(iNbLabelsRef-1 ));// The "-1" is to reject background label
    /*}
   else
   {
      if (m_)
      {
      }
      else
      {
      }
      po_fPPVL = 0;
      po_fSensL = std::numeric_limits<float>::infinity();
      bRes = false;
   }*/

    po_fF1 = 2 * (po_fPPVL * po_fSensL) / (po_fPPVL + po_fSensL);

    return bRes;
}

/**
@brief   Set the number of threads to use for the computation of ITK.
@param	[in] pi_iNbThreads Number of threads.
*/
void CAnalyzer::setNbThreads(int pi_iNbThreads)
{
    m_ThreadNb = static_cast<itk::ThreadIdType>(pi_iNbThreads);
}

/**
@brief  Compute a table of lesion overlap between 2 images.
@param	[out] po_iNbLabelsRef the number of labels found by connected components into reference image (1 by connected component + 1 for the background).
@param	[out] po_iNbLabelsTest the number of labels found by connected components into tested image (1 by connected component + 1 for the background).
@param	[out] po_ppiOverlapTab the table allocated by this method. Dimensions are po_iNbLabelsTest*po_iNbLabelsRef. The first cell [0][0] is the number of True-Negative voxels (background in both images), first line (excepte first cell) are False-Posite voxels, first colomne are False-Negative voxels and the rest of cels are True positive voxels. In Other words: first line is the back ground into Ref-image, first column is the back ground into Tested-image and the element [i][j] (with i,j > 0) correspond to the number of voxel belonging to the label i, that overlap with a voxel of the label j.
*/
void CAnalyzer::getOverlapTab(int&po_iNbLabelsRef, int&po_iNbLabelsTest, int* *&po_ppiOverlapTab)
{
    typedef itk::ImageRegionConstIterator<ImageType> imageConstIteratorType;
    typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> connectedComponentImageFilterType;
    typedef itk::RelabelComponentImageFilter<ImageType, ImageType> relabelComponentImageFilterType;
    typedef itk::ImageDuplicator<ImageType> imageDuplicatorFilterType;

    //Variable declaration
    connectedComponentImageFilterType::Pointer poLesionSeparatorFilter = connectedComponentImageFilterType::New();
    ImageType::Pointer poImageRefLesionsByLabels =  NULL;
    ImageType::Pointer poImageTestLesionsByLabels =  NULL;
    relabelComponentImageFilterType::Pointer poRelabelFilter = relabelComponentImageFilterType::New();
    int&iNbLabelsRef = po_iNbLabelsRef;
    int&iNbLabelsTest = po_iNbLabelsTest;
    int* *&ppiOverlapTab = po_ppiOverlapTab;

    //take into account the background label
    iNbLabelsRef = 1;
    iNbLabelsTest = 1;

    poLesionSeparatorFilter->SetNumberOfThreads(m_ThreadNb);

    //Create a label per connected component into ground truth
    poLesionSeparatorFilter->SetInput(m_imageRef);

    //Run filter
    poLesionSeparatorFilter->Update();

    // Compute Image Spacing, 4th dimension is not physical but temporal
    ImageType::SpacingType spacing = m_imageRef->GetSpacing();
    ImageType::SpacingValueType spacingTot = spacing[0];
    unsigned int imageDim = ImageType::ImageDimension;
    for (unsigned int i = 1; i < std::min(imageDim, (unsigned int)3); ++i)
    {
        spacingTot *= spacing[i];
    }

    // Compute minsize in voxels
    double minSizeInVoxelD = m_dfMinLesionVolumeDetection / spacingTot;
    double minSizeInVoxelD_floor = floor(minSizeInVoxelD);
    unsigned int minSizeInVoxel = static_cast<unsigned int>(minSizeInVoxelD_floor);
    minSizeInVoxel++; // to have strickly superior sizes

    poRelabelFilter->SetInput( poLesionSeparatorFilter->GetOutput() );
    poRelabelFilter->SetMinimumObjectSize(minSizeInVoxel);
    poRelabelFilter->SetNumberOfThreads(m_ThreadNb);
    poRelabelFilter->Update();

    iNbLabelsRef += poRelabelFilter->GetNumberOfObjects();
    //poImageRefLesionsByLabels =  poRelabelFilter->GetOutput();
    imageDuplicatorFilterType::Pointer poDuplicatorFilter = imageDuplicatorFilterType::New();
    poDuplicatorFilter->SetInputImage(poRelabelFilter->GetOutput());
    poDuplicatorFilter->Update();
    poImageRefLesionsByLabels = poDuplicatorFilter->GetOutput();


    //Create a label per connected component into image to evaluate
    poLesionSeparatorFilter->SetInput(m_imageTest);
    poLesionSeparatorFilter->Update();

    poRelabelFilter->SetInput( poLesionSeparatorFilter->GetOutput() );
    poRelabelFilter->SetMinimumObjectSize(minSizeInVoxel);
    poRelabelFilter->SetNumberOfThreads(m_ThreadNb);
    poRelabelFilter->Update();

    iNbLabelsTest += poRelabelFilter->GetNumberOfObjects();
    poImageTestLesionsByLabels =  poRelabelFilter->GetOutput();

    //Create and initialize table to store overlaps
    ppiOverlapTab = new int*[iNbLabelsRef];
    for(int i=0; i<iNbLabelsRef; ++i)
    {
        ppiOverlapTab[i] = new int[iNbLabelsTest];
        memset(ppiOverlapTab[i], 0, iNbLabelsTest*sizeof(int));
    }

    //Iterate on both image to fill the overlap tab
    imageConstIteratorType itRefLabels = imageConstIteratorType(poImageRefLesionsByLabels, poImageRefLesionsByLabels->GetLargestPossibleRegion());
    imageConstIteratorType itTestLabels = imageConstIteratorType(poImageTestLesionsByLabels, poImageTestLesionsByLabels->GetLargestPossibleRegion());

    while(!itRefLabels.IsAtEnd())
    {
        ImageType::PixelType voxelRefVal = itRefLabels.Value();
        ImageType::PixelType voxelTestVal = itTestLabels.Value();

        ppiOverlapTab[voxelRefVal][voxelTestVal]++;

        ++itRefLabels;
        ++itTestLabels;
    }

    return;
}

/**
   @brief  Getter of True positive lesions
   @return True positive lesions
*/
int CAnalyzer::getTruePositiveLesions(int pi_iNbLabelsRef, int pi_iNbLabelsTest, int * *pi_ppiOverlapTab)
{
    int iNbLesionsDetected = 0;

    double dfSensibility = 0;
    int *piTPRowSumTab = new int[pi_iNbLabelsRef]; //TODO verif si OK ou taille+1
    int *piColumnSumTab = new int[pi_iNbLabelsTest];

    //Init sum vectors
    memset(piTPRowSumTab, 0, pi_iNbLabelsRef*sizeof(int));
    memset(piColumnSumTab, 0, pi_iNbLabelsTest*sizeof(int));
    for (int i=0; i<pi_iNbLabelsRef; ++i) //TODO verif si < ou <=
    {
        for(int j=0; j<pi_iNbLabelsTest; ++j)
        {
            piColumnSumTab[j]+=pi_ppiOverlapTab[i][j];
            if (j>0)
            {
                //Iteration on each detection compared to a Ref label
                piTPRowSumTab[i] += pi_ppiOverlapTab[i][j];
            }
        }
    }

    //Iteration on each Ref label
    for(int i=1; i<pi_iNbLabelsRef; ++i)
    {
        int&iTP = piTPRowSumTab[i];  // True-Positive
        int&iFN = pi_ppiOverlapTab[i][0];  // False-Negative

        //Compute sensibility for one element of the ground-truth
        dfSensibility = ((double)iTP) / ((double)( iTP +  iFN ));

        //Check if sensibility is enough
        if (dfSensibility>m_dfDetectionThresholdAlpha)
        {
            //Call the second threshold test. Threshold on FalsePositive/DetectedVolume.
            if(falsePositiveRatioTester(i, pi_iNbLabelsTest, pi_iNbLabelsRef, pi_ppiOverlapTab, piTPRowSumTab, piColumnSumTab, m_dfDetectionThresholdBeta, m_dfDetectionThresholdGamma))
            {
                iNbLesionsDetected++;
            }
        }
    }

    return iNbLesionsDetected;
}

// Pair comparison functor
struct pair_decreasing_comparator
{
    bool operator() (const std::pair<int, int>&f, const std::pair<int, int>&s)
    {
        return (f.second > s.second);
    }
};


/**
@brief  Compute if a lesion is detected or not with beta and gamma thresholds.
@param	[out] pi_iLesionReference the number of labels found by connected components into reference image (1 by connected component + 1 for the background).
@param	[out] pi_iNbLabelsTest the number of labels found by connected components into tested image (1 by connected component + 1 for the background).
@param	[out] pi_iNbLabelsRef the number of the label of the tested lesion.
@param	[out] pi_ppiOverlapTab the mapping table.
@param	[out] pi_piTPRowSumTab the table of sum of Rows.
@param	[out] pi_piColumnSumTab the table of sum of Columns.
@param	[out] pi_dfBeta the beta threshold.
@param	[out] pi_dfGamma the beta threshold.
@return true if the lesion is detected, false in other cases.
*/
bool CAnalyzer::falsePositiveRatioTester(int pi_iLesionReference, int pi_iNbLabelsTest, int pi_iNbLabelsRef, int * *pi_ppiOverlapTab, int *pi_piTPRowSumTab, int *pi_piColumnSumTab, double pi_dfBeta, double pi_dfGamma)
{
    bool bRes = false;

    //////////////////////////////////////////////////////////////////////////
    // Locals variables declarations
    bool bExit = false;
    int k = 0;
    int &iSumOfTPForCurentRow = pi_piTPRowSumTab[pi_iLesionReference];
    double dfSumWeight = 0.0;
    double dfRatioOutsideInside = 0.0;   // This is the ratio of the outside part of a tested lesion on the total size od this tested lesion.
    std::vector<std::pair<int, int> > oSortedCollumVector;

    //////////////////////////////////////////////////////////////////////////
    // Construction of sorted vector for a TP column
    for(int l=1; l<pi_iNbLabelsTest; ++l)
    {
        int iTmpValOfTP = pi_ppiOverlapTab[pi_iLesionReference][l];
        oSortedCollumVector.push_back( std::pair<int, int>(l, iTmpValOfTP) );
    }
    std::sort(oSortedCollumVector.begin(), oSortedCollumVector.end(), pair_decreasing_comparator());

    //////////////////////////////////////////////////////////////////////////
    // Test in intersection size decreasing order that the regions overlapping the tested lesion are not too much outside of this lesion
    while(dfSumWeight<pi_dfGamma && k<pi_iNbLabelsRef && !bExit)
    {
        dfRatioOutsideInside = (double)(pi_ppiOverlapTab[0][oSortedCollumVector[k].first])/(double)(pi_piColumnSumTab[oSortedCollumVector[k].first]);
        bExit = dfRatioOutsideInside>pi_dfBeta;
        if (!bExit)
        {
            dfSumWeight += (double)(oSortedCollumVector[k].second)/(double)(iSumOfTPForCurentRow);
        }
        ++k;
    }

    //////////////////////////////////////////////////////////////////////////
    // Lesion is detected if a sufficient percentage (gamma) of the regions overlapping the tested lesion is not too much outside of this lesion
    bRes = !bExit; //Almost equivalent, except floating point imprecision, to bRes=dfSumWeight>=pi_dfGamma;

    return bRes;
}

/**
@brief  Compute a transposed table of lesion overlap between 2 images from the original table. It inverse reference image and tested image.
@param	[in]  pi_iNbLabelsRef the number of labels found by connected components into reference image (1 by connected component + 1 for the background).
@param	[in]  pi_iNbLabelsTest the number of labels found by connected components into tested image (1 by connected component + 1 for the background).
@param	[in]  pi_ppiOverlapTab the table allocated by this method. Dimensions are po_iNbLabelsTest*po_iNbLabelsRef.
@param	[out] po_rppiOverlapTabTransposed the transposed table allocated by this method. Dimensions are po_iNbLabelsTest*po_iNbLabelsRef.
*/
void CAnalyzer::transposer(int pi_iNbLabelsRef, int pi_iNbLabelsTest, int * *pi_ppiOverlapTab, int* *&po_rppiOverlapTabTransposed)
{
    po_rppiOverlapTabTransposed = new int*[pi_iNbLabelsTest+1];
    for(int i=0; i<pi_iNbLabelsTest; ++i)
    {
        po_rppiOverlapTabTransposed[i] = new int[pi_iNbLabelsRef+1];
    }

    for(int i=0; i<pi_iNbLabelsRef; ++i)
    {
        for(int j=0; j<pi_iNbLabelsTest; ++j)
        {
            po_rppiOverlapTabTransposed[j][i] = pi_ppiOverlapTab[i][j];
        }
    }
}

/**
@brief  Remove a table of lesion overlap between 2 images.
@param	[in] pi_ppiOverlapTab the table to remove.
@param	[in] pi_iNbLabelsRef the line dimension.
*/
void CAnalyzer::removeOverlapTab(int * *pi_ppiOverlapTab, int pi_iNbLabelsRef)
{
    if (pi_ppiOverlapTab)
    {
        for(int i=0; i<pi_iNbLabelsRef; ++i)
        {
            delete[] (pi_ppiOverlapTab[i]);
        }
        delete[] (pi_ppiOverlapTab);
    }
}

