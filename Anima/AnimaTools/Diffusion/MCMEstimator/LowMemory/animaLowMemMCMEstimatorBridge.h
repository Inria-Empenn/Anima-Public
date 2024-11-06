#pragma once

#include <animaImageDataSplitter.h>
#include <animaMCMEstimatorImageFilter.h>
#include <animaReadWriteFunctions.h>

namespace anima
{

class LowMemMCMEstimatorBridge
{
public:
    typedef anima::MCMEstimatorImageFilter <double, double> BaseFilterType;
    typedef BaseFilterType::InputImageType InputImageType;
    typedef BaseFilterType::OutputImageType OutputImageType;
    typedef BaseFilterType::OutputScalarImageType OutputScalarImageType;
    typedef BaseFilterType::MCMType MCMType;
    typedef BaseFilterType::SignalNoiseType SignalNoiseType;
    typedef BaseFilterType::MaximumLikelihoodEstimationMode MaximumLikelihoodEstimationMode;

    typedef anima::ImageDataSplitter < InputImageType > ImageSplitterType;
    typedef itk::Image <unsigned char,3> MaskImageType;
    typedef anima::ImageDataSplitter < MaskImageType > MoseImageSplitterType;

    LowMemMCMEstimatorBridge();
    ~LowMemMCMEstimatorBridge();

    std::string GetNameOfClass() {return "LowMemMCMEstimatorBridge";}

    void SetComputationMask(std::string &cMask);

    void SetInputMoseName(std::string &fileName);

    void SetDWIFileNames(std::string &fileList)
    {
        m_DWIImages->SetFileNames(fileList);
    }

    void SetOutputName(std::string &pref) {m_OutputName = pref;}

    void SetNbSplits(unsigned int nbSplits) {m_NbSplits = nbSplits;}
    void SetNumberOfThreads(unsigned int &nbT) {m_NumThreads = nbT;}

    void Update(int specificSplitToDo = -1, bool genOutputDescriptionData = false);

    void SetGradients(std::string grads) {m_Gradients = grads;}
    void SetBValues(std::string bvals) {m_BValues = bvals;}
    void SetBValueScale(bool val) {m_BValueScale = val;}

    void SetOutputAICName(std::string aicName) {m_OutputAICName = aicName;}
    void SetOutputB0Name(std::string b0Name) {m_OutputB0Name = b0Name;}
    void SetOutputSigmaName(std::string sigmaName) {m_OutputSigmaName = sigmaName;}
    void SetOutputMoseName(std::string moseName) {m_OutputMoseName = moseName;}

    void SetB0Threshold(double b0Thr) {m_B0Threshold = b0Thr;}
    void SetNumberOfFascicles(unsigned int nbFasc) {m_NumberOfFascicles = nbFasc;}
    void SetCompartmentType(unsigned int cType) {m_CompartmentType = cType;}

    void SetFreeWaterCompartment(bool fwComp) {m_FreeWaterCompartment = fwComp;}
    void SetStationaryWaterCompartment(bool swComp) {m_StationaryWaterCompartment = swComp;}
    void SetRestrictedWaterCompartment(bool irwComp) {m_RestrictedWaterCompartment = irwComp;}
    void SetStaniszCompartment(bool zComp) {m_StaniszCompartment = zComp;}

    void SetOptimizeFreeWaterDiffusivity(bool opt) {m_OptimizeFreeWaterDiffusivity = opt;}
    void SetOptimizeIRWDiffusivity(bool fix) {m_OptimizeIRWDiffusivity = fix;}
    void SetOptimizeStaniszDiffusivity(bool fix) {m_OptimizeStaniszDiffusivity = fix;}
    void SetOptimizeStaniszRadius(bool fix) {m_OptimizeStaniszRadius = fix;}

    void SetFixDiffusivity(bool fix) {m_FixDiffusivity = fix;}
    void SetFixKappa(bool fix) {m_FixKappa = fix;}
    void SetFixEAF(bool fix) {m_FixEAF = fix;}

    void SetAxialDiffusivityValue(double val) {m_AxialDiffusivityValue = val;}
    void SetRadialDiffusivity1Value(double val) {m_RadialDiffusivity1Value = val;}
    void SetRadialDiffusivity2Value(double val) {m_RadialDiffusivity2Value = val;}
    void SetIRWDiffusivityValue(double val) {m_IRWDiffusivityValue = val;}
    void SetStaniszDiffusivityValue(double val) {m_StaniszDiffusivityValue = val;}

    void SetCommonDiffusivities(bool common) {m_CommonDiffusivities = common;}
    void SetCommonKappa(bool common) {m_CommonKappa = common;}
    void SetCommonEAF(bool common) {m_CommonEAF = common;}

    void SetOptimizerType(std::string optType) {m_OptimizerType = optType;}
    void SetAbsCostChange(double num) {m_AbsCostChange = num;}
    void SetXTolerance(double num) {m_XTolerance = num;}
    void SetFTolerance(double num) {m_FTolerance = num;}
    void SetMaxEval(unsigned int num) {m_MaxEval = num;}
    void SetFindOptimalNumberOfCompartments(bool val) {m_FindOptimalNumberOfCompartments = val;}
    void SetMLEstimationStrategy(MaximumLikelihoodEstimationMode val) {m_MLEstimationStrategy = val;}

    void SetSmallDelta(double val) {m_SmallDelta = val;}
    void SetBigDelta(double val) {m_BigDelta = val;}

    void SetNoiseType(SignalNoiseType val) {m_NoiseType = val;}
    void SetNumberOfCoils(unsigned int val) {m_NumberOfCoils = val;}

    //Update progression of the process
    static void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
    {
        itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
        std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
    }

protected:
    void BuildAndWriteMCM(OutputImageType *tmpIm, std::string resName, OutputImageType::RegionType finalROI);

    template <typename TOutputType>
    void BuildAndWriteAdditional(TOutputType *tmpIm, std::string resName, typename TOutputType::RegionType finalROI)
    {
        typename TOutputType::RegionType tmpRegion = finalROI;
        for (unsigned int i = 0;i < TOutputType::GetImageDimension();++i)
            tmpRegion.SetIndex(i,0);

        typename TOutputType::Pointer tmpRes = TOutputType::New();
        tmpRes->Initialize();

        tmpRes->SetOrigin(m_ComputationMask->GetOrigin());
        tmpRes->SetRegions(tmpRegion);
        tmpRes->SetDirection(m_ComputationMask->GetDirection());
        tmpRes->SetSpacing(m_ComputationMask->GetSpacing());

        tmpRes->Allocate();

        itk::ImageRegionIterator <TOutputType> tmpImIt (tmpIm,finalROI);
        itk::ImageRegionIterator <TOutputType> tmpResIt (tmpRes,tmpRegion);

        while (!tmpImIt.IsAtEnd())
        {
            tmpResIt.Set(tmpImIt.Get());

            ++tmpImIt;
            ++tmpResIt;
        }

        anima::writeImage <TOutputType> (resName,tmpRes);
    }

private:
    std::string m_OutputName;

    unsigned int m_NbSplits;
    unsigned int m_NumThreads;

    ImageSplitterType *m_DWIImages;
    MoseImageSplitterType *m_InputMoseImage;
    MaskImageType::Pointer m_ComputationMask;

    std::string m_Gradients, m_BValues;
    bool m_BValueScale;
    std::string m_OutputAICName, m_OutputB0Name, m_OutputSigmaName, m_OutputMoseName;

    double m_B0Threshold;
    unsigned int m_NumberOfFascicles, m_CompartmentType;

    bool m_FreeWaterCompartment, m_StationaryWaterCompartment, m_RestrictedWaterCompartment, m_StaniszCompartment;
    MaximumLikelihoodEstimationMode m_MLEstimationStrategy;
    bool m_FindOptimalNumberOfCompartments;

    bool m_FixDiffusivity, m_OptimizeFreeWaterDiffusivity, m_OptimizeIRWDiffusivity, m_FixKappa, m_FixEAF;
    bool m_OptimizeStaniszDiffusivity, m_OptimizeStaniszRadius;
    bool m_CommonDiffusivities, m_CommonKappa, m_CommonEAF;

    double m_AxialDiffusivityValue;
    double m_RadialDiffusivity1Value, m_RadialDiffusivity2Value;
    double m_IRWDiffusivityValue;
    double m_StaniszDiffusivityValue;

    std::string m_OptimizerType;
    double m_AbsCostChange;
    double m_XTolerance, m_FTolerance;
    unsigned int m_MaxEval;

    double m_SmallDelta;
    double m_BigDelta;

    SignalNoiseType m_NoiseType;
    unsigned int m_NumberOfCoils;
};

} // end namespace anima
