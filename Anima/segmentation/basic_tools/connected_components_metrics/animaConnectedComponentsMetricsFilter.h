#pragma once

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include "animaConnectedComponentsVolumeFilter.h"
#include <animaReadWriteFunctions.h>

namespace anima
{
/** @brief Compute connected components metrics.
 *
 *
 */
template <typename ImageType, typename OutputType = ImageType>
class ConnectedComponentsMetricsFilter :
public itk::ImageToImageFilter< ImageType , OutputType >
{
public:
    /** Standard class typedefs. */
    typedef ConnectedComponentsMetricsFilter Self;
    typedef itk::ImageToImageFilter< ImageType , OutputType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(ConnectedComponentsMetricsFilter, ImageToImageFilter);

    /** Image typedef support */

    /**  Type of image. */
    typedef typename ImageType::ConstPointer ImageConstPointerType;
    typedef typename itk::ImageRegionConstIterator< ImageType > ImageConstIteratorType;

    typedef typename OutputType::Pointer OutputImagePointerType;
    typedef typename OutputType::PixelType OutputPixelType;
    typedef typename itk::ImageRegionIterator< OutputType > OutputIteratorType;

    typedef anima::ConnectedComponentVolumeFilter<ImageType> ConnectedComponentVolumeFilterType;

    typedef typename ImageType::SpacingValueType spacingValueType;

    void WriteOutputs();

    void SetInputReference(const ImageType* image);
    void SetInputTest(const ImageType* image);

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
        m_Tol = tol;
    }

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(MinSizeMM3, double);
    itkGetMacro(MinSizeMM3, double);

    itkSetMacro(Beta, double);
    itkGetMacro(Beta, double);

    itkSetMacro(Alpha, double);
    itkGetMacro(Alpha, double);

    itkSetMacro(FullyConnected, bool);
    itkGetMacro(FullyConnected, bool);

    itkSetMacro(Verbose, bool);
    itkGetMacro(Verbose, bool);

    itkSetMacro(OverlapDetectionType, bool);
    itkGetMacro(OverlapDetectionType, bool);

    itkSetMacro(Dimension, unsigned int);
    itkGetMacro(Dimension, unsigned int);

    void SetOutputFilenameCCRef(std::string outputFilename) {m_OutputCCReference_filename=outputFilename;}
    void SetOutputFilenameCCTest(std::string outputFilename) {m_OutputCCTest_filename=outputFilename;}

    void SetOutputFilenameDiffVoxelWise(std::string outputFilename) {m_OutputDiffVoxelWise_filename=outputFilename;}
    void SetOutputFilenameDiffTest(std::string outputFilename) {m_OutputDiffTest_filename=outputFilename;}
    void SetOutputFilenameDiffRef(std::string outputFilename) {m_OutputDiffRef_filename=outputFilename;}

    void SetOutputFilenameEvolutionTest(std::string outputFilename) {m_OutputEvolutionTest_filename=outputFilename;}
    void SetOutputFilenameEvolutionRef(std::string outputFilename) {m_OutputEvolutionRef_filename=outputFilename;}

    void SetOutputFilenameTextVolTest(std::string outputFilename) {m_volume_test_text_filename=outputFilename;}
    void SetOutputFilenameTextVolRef(std::string outputFilename) {m_volume_ref_text_filename=outputFilename;}
    void SetOutputFilenameTextMetrics(std::string outputFilename) {m_metrics_text_filename=outputFilename;}

    itkGetMacro(OriginalNumberOfObjectsTest, unsigned int);
    itkGetMacro(NumberOfObjectsTest, unsigned int);
    itkGetMacro(TestTotalVolume, double);

    itkGetMacro(OriginalNumberOfObjectsRef, unsigned int);
    itkGetMacro(NumberOfObjectsRef, unsigned int);
    itkGetMacro(ReferenceTotalVolume, double);

    itkGetMacro(nb_true_positive_test, unsigned int);
    itkGetMacro(nb_true_positive_ref, unsigned int);
    itkGetMacro(nb_false_positive, unsigned int);
    itkGetMacro(nb_false_negative, unsigned int);

    itkGetMacro(nb_grow_ref, unsigned int);
    itkGetMacro(nb_shrink_ref, unsigned int);
    itkGetMacro(nb_same_ref, unsigned int);

    itkGetMacro(nb_grow_test, unsigned int);
    itkGetMacro(nb_shrink_test, unsigned int);
    itkGetMacro(nb_same_test, unsigned int);

    itkGetMacro(VolumeDifference, double);
    itkGetMacro(VolumeDifferenceRatio, double);

    itkGetMacro(FalsePositiveRate, double);
    itkGetMacro(FalseNegativeRate, double);
    itkGetMacro(RefTruePositiveRate, double);
    itkGetMacro(TestTruePositiveRate, double);

    itkGetMacro(RefGrowRate, double);
    itkGetMacro(RefShrinkRate, double);
    itkGetMacro(RefSameRate, double);

    itkGetMacro(TestGrowRate, double);
    itkGetMacro(TestShrinkRate, double);
    itkGetMacro(TestSameRate, double);

    std::vector<unsigned int> GetGrownTest(){return m_vect_grow_test;}
    std::vector<unsigned int> GetShrinkTest(){return m_vect_shrink_test;}
    std::vector<unsigned int> GetSameTest(){return m_vect_same_test;}

    std::vector<unsigned int> GetGrownRef(){return m_vect_grow_ref;}
    std::vector<unsigned int> GetShrinkRef(){return m_vect_shrink_ref;}
    std::vector<unsigned int> GetSameRef(){return m_vect_same_ref;}

    std::vector<unsigned int> GetTestTP(){return m_vect_true_positive_test;}
    std::vector<unsigned int> GetRefTP(){return m_vect_true_positive_ref;}

    std::vector<double> GetVectVolumeTest(){return m_vect_volume_test;}
    std::vector<double> GetVectVolumeReference(){return m_vect_volume_ref;}

    OutputType* GetOutputCCReference();
    OutputType* GetOutputCCTest();
    OutputType* GetOutputDiffVoxelWise();
    OutputType* GetOutputDiffTest();
    OutputType* GetOutputDiffRef();

    OutputType* GetOutputEvolutionRef();
    OutputType* GetOutputEvolutionTest();

protected:

    ConnectedComponentsMetricsFilter()
    {
        this->SetNumberOfRequiredOutputs(7);
        this->SetNumberOfRequiredInputs(2);

        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );
        this->SetNthOutput( 2, this->MakeOutput(2) );
        this->SetNthOutput( 3, this->MakeOutput(3) );
        this->SetNthOutput( 4, this->MakeOutput(4) );
        this->SetNthOutput( 5, this->MakeOutput(5) );
        this->SetNthOutput( 6, this->MakeOutput(6) );



        /**
         * Input parameters
         * */
        m_OverlapDetectionType = false;
        m_MinSizeMM3 = 0;
        m_FullyConnected = false;
        m_Verbose = false;
        m_Tol = 0.0001;
        m_Beta = 1;
        m_Alpha = 1;
        m_SpacingTot = 1;
        m_Dimension = 2;
        this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());

        /**
         * Images features
         * */
        m_OriginalNumberOfObjectsTest = 0;
        m_NumberOfObjectsTest = 0;
        m_TestTotalVolume = 0;

        m_OriginalNumberOfObjectsRef = 0;
        m_NumberOfObjectsRef = 0;
        m_ReferenceTotalVolume = 0;


        /**
         * Detection
         * */
        m_nb_true_positive_test = 0;
        m_nb_true_positive_ref = 0;
        m_nb_false_positive = 0;
        m_nb_false_negative = 0;

        /**
         * Evolution
         * */
        m_nb_grow_ref = 0;
        m_nb_shrink_ref = 0;
        m_nb_same_ref = 0;

        m_nb_grow_test = 0;
        m_nb_shrink_test = 0;
        m_nb_same_test = 0;

        /**
         * Metrics
         * */
        m_VolumeDifference = 0;
        m_VolumeDifferenceRatio = std::numeric_limits<double>::quiet_NaN();

        m_FalsePositiveRate = std::numeric_limits<double>::quiet_NaN();
        m_FalseNegativeRate = std::numeric_limits<double>::quiet_NaN();
        m_RefTruePositiveRate = std::numeric_limits<double>::quiet_NaN();
        m_TestTruePositiveRate = std::numeric_limits<double>::quiet_NaN();

        m_RefGrowRate = std::numeric_limits<double>::quiet_NaN();
        m_RefShrinkRate = std::numeric_limits<double>::quiet_NaN();
        m_RefSameRate = std::numeric_limits<double>::quiet_NaN();

        m_TestGrowRate = std::numeric_limits<double>::quiet_NaN();
        m_TestShrinkRate = std::numeric_limits<double>::quiet_NaN();
        m_TestSameRate = std::numeric_limits<double>::quiet_NaN();
    }

    virtual ~ConnectedComponentsMetricsFilter()
    {
    }

    itk::DataObject::Pointer MakeOutput(unsigned int idx);

    typename ImageType::ConstPointer GetInputReference();
    typename ImageType::ConstPointer GetInputTest();

    void PrintInfo();
    void PrintMetricsFile();
    void ComputeMetrics();
    void ComputeIntersectionBetweenCC();
    void ComputeEvolution();
    void ComputeDetection();
    void GenerateData();
    void PrintImagesOutput();

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

private:
    ConnectedComponentsMetricsFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /**
     * Input parameters
     * */
    bool m_OverlapDetectionType;
    double m_MinSizeMM3;
    bool m_FullyConnected;
    unsigned int m_Dimension;
    bool m_Verbose;
    double m_Tol;
    double m_Beta;
    double m_Alpha;
    spacingValueType m_SpacingTot;

    /**
     * Output filenames
     * */

    // image file
    std::string m_OutputCCReference_filename;
    std::string m_OutputCCTest_filename;
    std::string m_OutputDiffVoxelWise_filename;
    std::string m_OutputDiffRef_filename;
    std::string m_OutputDiffTest_filename;
    std::string m_OutputEvolutionTest_filename;
    std::string m_OutputEvolutionRef_filename;

    // information text file
    std::string m_volume_test_text_filename;
    std::string m_volume_ref_text_filename;
    std::string m_metrics_text_filename;

    /**
     * Images volumes
     * */
    unsigned int m_OriginalNumberOfObjectsTest;                             /*!< Original number of objects in test image */
    unsigned int m_NumberOfObjectsTest;                                     /*!< Number of objects in test image (after removing small ones) */
    std::vector<double> m_vect_volume_test;                                 /*!< Volume of each object */
    double m_TestTotalVolume;                                               /*!< Total volume of test image */

    unsigned int m_OriginalNumberOfObjectsRef;                              /*!< Original number of objects in reference image */
    unsigned int m_NumberOfObjectsRef;                                      /*!< Number of objects in reference image (after removing small ones) */
    std::vector<double> m_vect_volume_ref;                                  /*!< Volume of each object */
    double m_ReferenceTotalVolume;                                          /*!< Total volume of reference image */

    /**
     * Interaction between connected components
     * */
    std::vector<std::vector<unsigned short int> > m_IntersectionVoxels;
    std::vector<std::vector<double> > m_IntersectionMM3;
    std::vector<std::vector<unsigned short int> > m_RefInclusDansTest;
    std::vector<std::vector<unsigned short int> > m_TestInclusDansRef;

    /**
     * Detection
     * */
    std::vector<unsigned short int> m_vect_true_positive_test;              /*!< List of true positive objects in test image */
    std::vector<unsigned short int> m_vect_true_positive_ref;               /*!< List of true positive objects in reference image */

    unsigned int m_nb_true_positive_test;                                   /*!< Number of true positive objects in test image */
    unsigned int m_nb_true_positive_ref;                                    /*!< Number of true positive objects in reference image */
    unsigned int m_nb_false_positive;                                       /*!< Number of false positive objects (in test image) */
    unsigned int m_nb_false_negative;                                       /*!< Number of false negative objects (in reference image) */

    /**
     * Evolution
     * */

    unsigned int m_nb_grow_ref;                                              /*!< Number of objects that will grow (in reference image) */
    unsigned int m_nb_shrink_ref;                                            /*!< Number of objects that will shrink (in reference image) */
    unsigned int m_nb_same_ref;                                              /*!< Number of objects that stay the same (in reference image) */

    unsigned int m_nb_grow_test;                                             /*!< Number of objects have grown (in test image) */
    unsigned int m_nb_shrink_test;                                           /*!< Number of objects have shrinked (in test image) */
    unsigned int m_nb_same_test;                                             /*!< Number of objects that stay the same (in test image) */

    std::vector<unsigned short int> m_vect_grow_test;                        /*!< List of objects that will grow (in reference image) */
    std::vector<unsigned short int> m_vect_shrink_test;                      /*!< List of objects that will shrink (in reference image) */
    std::vector<unsigned short int> m_vect_same_test;                        /*!< List of objects that stay the same (in reference image) */

    std::vector<unsigned short int> m_vect_grow_ref;                         /*!< List of objects have grown (in test image) */
    std::vector<unsigned short int> m_vect_shrink_ref;                       /*!< List of objects have shrinked (in test image) */
    std::vector<unsigned short int> m_vect_same_ref;                         /*!< List of objects that stay the same (in test image) */


    /**
     * Metrics
     * */

    double m_VolumeDifference;                                                /*!< Volume Difference –  absolute difference in volumes */
    double m_VolumeDifferenceRatio;                                           /*!< Volume Difference Ratio –  absolute difference in volumes divided by the true volume */

    double m_FalsePositiveRate;                                               /*!< Object wise false positive rate */
    double m_FalseNegativeRate;                                               /*!< Object wise false negative rate */
    double m_RefTruePositiveRate;                                             /*!< Object wise true positive rate in reference image (= recall) */
    double m_TestTruePositiveRate;                                            /*!< Object wise true positive rate in test image (= precision) */

    double m_RefGrowRate;                                                     /*!< Rate of growing objects in reference image */
    double m_RefShrinkRate;                                                   /*!< Rate of shrinking objects in reference image */
    double m_RefSameRate;                                                     /*!< Rate of same objects in reference image */

    double m_TestGrowRate;                                                    /*!< Rate of growing objects in test image */
    double m_TestShrinkRate;                                                  /*!< Rate of shrinking objects in test image */
    double m_TestSameRate;                                                    /*!< Rate of same objects in test image */

};
} // end of namespace anima

#include "animaConnectedComponentsMetricsFilter.hxx"

