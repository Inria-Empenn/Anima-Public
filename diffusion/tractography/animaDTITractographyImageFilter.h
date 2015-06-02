#pragma once

#include <itkVectorImage.h>
#include <itkImage.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkProcessObject.h>

#include "AnimaTractographyExport.h"

#include <vector>

namespace anima
{

class ANIMATRACTOGRAPHY_EXPORT dtiTractographyImageFilter : public itk::ProcessObject
{
public:
    /** SmartPointer typedef support  */
    typedef dtiTractographyImageFilter Self;
    typedef itk::ProcessObject Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self);

    itkTypeMacro(dtiTractographyImageFilter,itk::ProcessObject);

    typedef enum {
        Forward,
        Backward,
        Both
    } FiberProgressType;
    
    typedef itk::VectorImage <float, 3> LogDTIImageType;
    typedef LogDTIImageType::Pointer LogDTIImagePointer;
    typedef LogDTIImageType::PixelType VectorType;
    typedef itk::LinearInterpolateImageFunction <LogDTIImageType> DTIInterpolatorType;
    typedef DTIInterpolatorType::Pointer DTIInterpolatorPointer;
    typedef itk::ContinuousIndex <double, 3> ContinuousIndexType;
    
    typedef itk::Image <unsigned short, 3> MaskImageType;
    typedef MaskImageType::Pointer MaskImagePointer;
    typedef MaskImageType::PointType PointType;
    typedef MaskImageType::IndexType IndexType;
    
    typedef std::vector <PointType> FiberType;
    typedef std::vector <std::pair < FiberProgressType, FiberType > > FiberProcessVectorType;
    
    typedef struct {
        dtiTractographyImageFilter *trackerPtr;
        std::vector < std::vector <FiberType> > resultFibersFromThreads;
    } trackerArguments;
    
    void SetInputImage(LogDTIImageType *input) {m_InputImage = input;}
    void SetTrackingMask(MaskImageType *mask) {m_MaskImage = mask;}
    void SetForbiddenMask(MaskImageType *mask) {m_ForbiddenMaskImage = mask;}
    void SetCutMask(MaskImageType *mask) {m_CutMaskImage = mask;}
    
    void SetNumberOfThreads(unsigned int numThr) {m_NumberOfThreads = numThr;}
    void SetNumberOfFibersPerPixel(unsigned int num) {m_NumberOfFibersPerPixel = num;}
    
    void SetStepProgression(double num) {m_StepProgression = num;}
    void SetStopFAThreshold(double num) {m_StopFAThreshold = num;}
    
    void SetMaxFiberAngle(double num) {m_MaxFiberAngle = num;}
    void SetMinLengthFiber(double num) {m_MinLengthFiber = num;}
    void SetMaxLengthFiber(double num) {m_MaxLengthFiber = num;}
    
    void Update();
    
    void SetComputeLocalColors(bool flag) {m_ComputeLocalColors = flag;}
    void createVTKOutput(std::vector < std::vector <PointType> > &filteredFibers);
    vtkPolyData *GetOutput() {return m_Output;}
    
protected:
    dtiTractographyImageFilter();
    virtual ~dtiTractographyImageFilter();

    static ITK_THREAD_RETURN_TYPE ThreadTracker(void *arg);
    void ThreadTrack(unsigned int numThread, std::vector <FiberType> &resultFibers);
    
    FiberProcessVectorType ComputeFiber(FiberType &fiber, DTIInterpolatorType * dtiInterpolator,
                                        FiberProgressType ways);
    
    void PrepareTractography();
    std::vector < FiberType > FilterOutputFibers(std::vector < FiberType > &fibers);
    
    double GetFA(VectorType &dtiLogValue);
    std::vector <double> GetDTIPrincipalDirection(VectorType &dtiLogValue, bool is2d);
    
private:
    unsigned int m_NumberOfThreads;
    unsigned int m_NumberOfFibersPerPixel;
    
    double m_StepProgression;
    double m_StopFAThreshold;
    double m_MaxFiberAngle;
    
    double m_MinLengthFiber;
    double m_MaxLengthFiber;
    
    LogDTIImagePointer m_InputImage;
    MaskImagePointer m_MaskImage, m_ForbiddenMaskImage, m_CutMaskImage;
    
    FiberProcessVectorType m_PointsToProcess;
    std::vector <unsigned int> m_FilteringValues;
    
    bool m_ComputeLocalColors;
    vtkSmartPointer<vtkPolyData> m_Output;
};

} // end of namespace anima
