#include <animaMultiCompartmentModelCreator.h>

#include <animaFreeWaterCompartment.h>
#include <animaIsotropicRestrictedWaterCompartment.h>
#include <animaNODDICompartment.h>
#include <animaStationaryWaterCompartment.h>
#include <animaStaniszCompartment.h>
#include <animaStickCompartment.h>
#include <animaTensorCompartment.h>
#include <animaZeppelinCompartment.h>

namespace anima
{

MultiCompartmentModelCreator::MultiCompartmentModelCreator()
{
    m_CompartmentType = Tensor;
    m_ModelWithFreeWaterComponent = true;
    m_ModelWithStationaryWaterComponent = false;
    m_ModelWithRestrictedWaterComponent = false;
    m_ModelWithStaniszComponent = false;

    m_FreeWaterProportionFixedValue = 0;
    m_StationaryWaterProportionFixedValue = 0;
    m_RestrictedWaterProportionFixedValue = 0;
    m_StaniszProportionFixedValue = 0;

    m_NumberOfCompartments = 1;

    m_UseFixedWeights = false;
    m_UseConstrainedDiffusivity = false;
    m_UseConstrainedFreeWaterDiffusivity = true;
    m_UseConstrainedIRWDiffusivity = false;
    m_UseConstrainedOrientationConcentration = false;
    m_UseConstrainedExtraAxonalFraction = false;
    m_UseBoundedOptimization = false;

    m_UseCommonDiffusivities = false;
    m_UseCommonConcentrations = false;
    m_UseCommonExtraAxonalFractions = false;

    m_AxialDiffusivity = 1.71e-3;
    m_FreeWaterDiffusivity = 3.0e-3;
    m_RadialDiffusivity1 = 1.5e-4;
    m_RadialDiffusivity2 = 1.5e-4;
    m_ExtraAxonalFraction = 0.1;
    m_OrientationConcentration = 10.0;

    m_UserDefinedConcentrationBounds = false;
    m_ConcentrationLowerBound = 0;
    m_ConcentrationUpperBound = 0;
}

void MultiCompartmentModelCreator::SetConcentrationBounds(double lowerBound, double upperBound)
{
    m_UserDefinedConcentrationBounds = true;
    m_ConcentrationLowerBound = lowerBound;
    m_ConcentrationUpperBound = upperBound;
    m_OrientationConcentration = (lowerBound + upperBound) / 2.0;
}

MultiCompartmentModelCreator::MCMPointer MultiCompartmentModelCreator::GetNewMultiCompartmentModel()
{
    MCMPointer outputMCM = MCMType::New();
    outputMCM->SetOptimizeWeights(!m_UseFixedWeights);
    outputMCM->SetUseBoundedWeightsOptimization(m_UseBoundedOptimization);
    outputMCM->SetCommonDiffusivityParameters(m_UseCommonDiffusivities);
    outputMCM->SetCommonConcentrationParameters(m_UseCommonConcentrations);
    outputMCM->SetCommonExtraAxonalFractionParameters(m_UseCommonExtraAxonalFractions);

    double sumIsotropicWeights = 0;
    if (m_ModelWithFreeWaterComponent)
        sumIsotropicWeights += m_FreeWaterProportionFixedValue;

    if (m_ModelWithStationaryWaterComponent)
        sumIsotropicWeights += m_StationaryWaterProportionFixedValue;

    if (m_ModelWithRestrictedWaterComponent)
        sumIsotropicWeights += m_RestrictedWaterProportionFixedValue;

    if (m_ModelWithStaniszComponent)
        sumIsotropicWeights += m_StaniszProportionFixedValue;

    double compartmentWeight = (1.0 - sumIsotropicWeights) / m_NumberOfCompartments;

    if (m_ModelWithFreeWaterComponent)
    {
        typedef anima::FreeWaterCompartment FreeWaterType;
        FreeWaterType::Pointer fwComp = FreeWaterType::New();
        fwComp->SetEstimateAxialDiffusivity(!m_UseConstrainedFreeWaterDiffusivity);
        fwComp->SetAxialDiffusivity(m_FreeWaterDiffusivity);
        fwComp->SetUseBoundedOptimization(m_UseBoundedOptimization);

        if (m_NumberOfCompartments != 0)
            outputMCM->AddCompartment(m_FreeWaterProportionFixedValue,fwComp);
        else
            outputMCM->AddCompartment(m_FreeWaterProportionFixedValue / sumIsotropicWeights,fwComp);
    }

    if (m_ModelWithStationaryWaterComponent)
    {
        typedef anima::StationaryWaterCompartment SWType;
        SWType::Pointer swComp = SWType::New();

        if (m_NumberOfCompartments != 0)
            outputMCM->AddCompartment(m_StationaryWaterProportionFixedValue,swComp);
        else
            outputMCM->AddCompartment(m_StationaryWaterProportionFixedValue / sumIsotropicWeights,swComp);
    }

    if (m_ModelWithRestrictedWaterComponent)
    {
        typedef anima::IsotropicRestrictedWaterCompartment IRWType;
        IRWType::Pointer restComp = IRWType::New();
        restComp->SetEstimateAxialDiffusivity(!m_UseConstrainedIRWDiffusivity);
        restComp->SetUseBoundedOptimization(m_UseBoundedOptimization);

        if (m_NumberOfCompartments != 0)
            outputMCM->AddCompartment(m_RestrictedWaterProportionFixedValue,restComp);
        else
            outputMCM->AddCompartment(m_RestrictedWaterProportionFixedValue / sumIsotropicWeights,restComp);
    }

    if (m_ModelWithStaniszComponent)
    {
        typedef anima::StaniszCompartment StaniszType;
        StaniszType::Pointer restComp = StaniszType::New();
        restComp->SetEstimateAxialDiffusivity(!m_UseConstrainedDiffusivity);
        restComp->SetUseBoundedOptimization(m_UseBoundedOptimization);

        if (m_NumberOfCompartments != 0)
            outputMCM->AddCompartment(m_StaniszProportionFixedValue,restComp);
        else
            outputMCM->AddCompartment(m_StaniszProportionFixedValue / sumIsotropicWeights,restComp);
    }

    for (unsigned int i = 0;i < m_NumberOfCompartments;++i)
    {
        anima::BaseCompartment::Pointer tmpPointer;
        bool applyCommonConstraints = (i > 0);

        switch (m_CompartmentType)
        {
            case Stick:
                this->CreateStickCompartment(tmpPointer,applyCommonConstraints);
                break;

            case Zeppelin:
                this->CreateZeppelinCompartment(tmpPointer,applyCommonConstraints);
                break;

            case Tensor:
                this->CreateTensorCompartment(tmpPointer,applyCommonConstraints);
                break;
                
            case NODDI:
                this->CreateNODDICompartment(tmpPointer,applyCommonConstraints);
                break;

            case DDI:
                this->CreateDDICompartment(tmpPointer,applyCommonConstraints);
                break;

            default:
                throw itk::ExceptionObject(__FILE__, __LINE__,"Creation of multiple free water compartment model not handled",ITK_LOCATION);
                break;
        }

        tmpPointer->SetUseBoundedOptimization(m_UseBoundedOptimization);
        // Kind of ugly but required for optimization, otherwise initialization from simplified models may fail
        tmpPointer->SetOrientationConcentration(m_OrientationConcentration);

        outputMCM->AddCompartment(compartmentWeight,tmpPointer);
    }

    return outputMCM;
}

void MultiCompartmentModelCreator::CreateStickCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
{
    typedef anima::StickCompartment StickType;

    StickType::Pointer stickComp = StickType::New();
    stickComp->SetEstimateAxialDiffusivity(!m_UseConstrainedDiffusivity);

    stickComp->SetAxialDiffusivity(m_AxialDiffusivity);
    stickComp->SetRadialDiffusivity1((m_RadialDiffusivity1 + m_RadialDiffusivity2) / 2.0);
    
    if (applyConstraints)
    {
        if (m_UseCommonDiffusivities)
            stickComp->SetEstimateAxialDiffusivity(false);
    }

    compartmentPointer = stickComp;
}

void MultiCompartmentModelCreator::CreateZeppelinCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
{
    typedef anima::ZeppelinCompartment ZeppelinType;

    ZeppelinType::Pointer zepComp = ZeppelinType::New();
    zepComp->SetEstimateDiffusivities(!m_UseConstrainedDiffusivity);

    zepComp->SetAxialDiffusivity(m_AxialDiffusivity);
    zepComp->SetRadialDiffusivity1((m_RadialDiffusivity1 + m_RadialDiffusivity2) / 2.0);

    if (applyConstraints)
    {
        if (m_UseCommonDiffusivities)
            zepComp->SetEstimateDiffusivities(false);
    }

    compartmentPointer = zepComp;
}

void MultiCompartmentModelCreator::CreateTensorCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
{
    typedef anima::TensorCompartment TensorType;

    TensorType::Pointer tensComp = TensorType::New();
    tensComp->SetEstimateDiffusivities(!m_UseConstrainedDiffusivity);

    tensComp->SetAxialDiffusivity(m_AxialDiffusivity);
    tensComp->SetRadialDiffusivity1(m_RadialDiffusivity1);
    tensComp->SetRadialDiffusivity2(m_RadialDiffusivity2);

    if (applyConstraints)
    {
        if (m_UseCommonDiffusivities)
            tensComp->SetEstimateDiffusivities(false);
    }

    compartmentPointer = tensComp;
}

void MultiCompartmentModelCreator::CreateNODDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
{
    typedef anima::NODDICompartment NODDIType;
    
    NODDIType::Pointer noddiComp = NODDIType::New();
    noddiComp->SetEstimateAxialDiffusivity(!m_UseConstrainedDiffusivity);
    
    noddiComp->SetOrientationConcentration(m_OrientationConcentration);
    noddiComp->SetExtraAxonalFraction(m_ExtraAxonalFraction);
    noddiComp->SetAxialDiffusivity(m_AxialDiffusivity);
    
    if (applyConstraints)
    {
        if (m_UseCommonDiffusivities)
            noddiComp->SetEstimateAxialDiffusivity(false);

        if (this->GetUseCommonConcentrations())
            noddiComp->SetEstimateOrientationConcentration(false);

        if (this->GetUseCommonExtraAxonalFractions())
            noddiComp->SetEstimateExtraAxonalFraction(false);
    }
    
    compartmentPointer = noddiComp;
}

void MultiCompartmentModelCreator::CreateDDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
{
    std::string error("DDI model not implemented in the public version of ANIMA");
    throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
}

} // end namespace anima
