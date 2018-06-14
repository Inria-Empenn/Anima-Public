#pragma once

#include <animaMultiCompartmentModel.h>
#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>

namespace anima
{

//! Really this class is some simplified factory that creates the MCM that it knows
// One could think of making it a singleton, but parameter are somewhat in the way
class ANIMAMCM_EXPORT MultiCompartmentModelCreator
{
public:
    MultiCompartmentModelCreator();
    virtual ~MultiCompartmentModelCreator() {}

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Pointer MCMPointer;

    typedef anima::BaseCompartment BaseCompartmentType;
    typedef BaseCompartmentType::Pointer BaseCompartmentPointer;
    typedef anima::DiffusionModelCompartmentType CompartmentType;

    // Isotropic compartments (one free water, one stationary water, one restricted water
    void SetModelWithFreeWaterComponent(bool arg) {m_ModelWithFreeWaterComponent = arg;}
    void SetModelWithStationaryWaterComponent(bool arg) {m_ModelWithStationaryWaterComponent = arg;}
    void SetModelWithRestrictedWaterComponent(bool arg) {m_ModelWithRestrictedWaterComponent = arg;}
    void SetFreeWaterProportionFixedValue(double arg) {m_FreeWaterProportionFixedValue = arg;}
    void SetStationaryWaterProportionFixedValue(double arg) {m_StationaryWaterProportionFixedValue = arg;}
    void SetRestrictedWaterProportionFixedValue(double arg) {m_RestrictedWaterProportionFixedValue = arg;}

    void SetCompartmentType(CompartmentType arg) {m_CompartmentType = arg;}
    void SetNumberOfCompartments(unsigned int num) {m_NumberOfCompartments = num;}
    void SetUseFixedWeights(bool arg) {m_UseFixedWeights = arg;}

    void SetUseConstrainedDiffusivity(bool arg) {m_UseConstrainedDiffusivity = arg;}
    void SetUseConstrainedOrientationConcentration(bool arg) {m_UseConstrainedOrientationConcentration = arg;}
    void SetUseConstrainedExtraAxonalFraction(bool arg) {m_UseConstrainedExtraAxonalFraction = arg;}
    void SetUseConstrainedFreeWaterDiffusivity(bool arg) {m_UseConstrainedFreeWaterDiffusivity = arg;}
    void SetUseConstrainedIRWDiffusivity(bool arg) {m_UseConstrainedIRWDiffusivity = arg;}
    void SetUseBoundedOptimization(bool arg) {m_UseBoundedOptimization = arg;}

    bool GetUseConstrainedDiffusivity() {return m_UseConstrainedDiffusivity;}
    bool GetUseConstrainedOrientationConcentration() {return m_UseConstrainedOrientationConcentration;}
    bool GetUseConstrainedExtraAxonalFraction() {return m_UseConstrainedExtraAxonalFraction;}

    void SetUseCommonDiffusivities(bool arg) {m_UseCommonDiffusivities = arg;}
    void SetUseCommonConcentrations(bool arg) {m_UseCommonConcentrations = arg;}
    void SetUseCommonExtraAxonalFractions(bool arg) {m_UseCommonExtraAxonalFractions = arg;}

    bool GetUseCommonDiffusivities() {return m_UseCommonDiffusivities;}
    bool GetUseCommonConcentrations() {return m_UseCommonConcentrations;}
    bool GetUseCommonExtraAxonalFractions() {return m_UseCommonExtraAxonalFractions;}

    void SetFreeWaterDiffusivityValue(double arg) {m_FreeWaterDiffusivity = arg;}
    void SetAxialDiffusivityValue(double arg) {m_AxialDiffusivity = arg;}
    void SetRadialDiffusivity1Value(double arg) {m_RadialDiffusivity1 = arg;}
    void SetRadialDiffusivity2Value(double arg) {m_RadialDiffusivity2 = arg;}

    void SetConcentrationBounds(double lowerBound, double upperBound);
    bool GetUserDefinedConcentrationBounds() {return m_UserDefinedConcentrationBounds;}

    double GetOrientationConcentration() {return m_OrientationConcentration;}
    double GetAxialDiffusivity() {return m_AxialDiffusivity;}
    double GetRadialDiffusivity1() {return m_RadialDiffusivity1;}
    double GetRadialDiffusivity2() {return m_RadialDiffusivity2;}
    double GetExtraAxonalFraction() {return m_ExtraAxonalFraction;}
    double GetConcentrationLowerBound() {return m_ConcentrationLowerBound;}
    double GetConcentrationUpperBound() {return m_ConcentrationUpperBound;}

    MCMPointer GetNewMultiCompartmentModel();

private:
    void CreateStickCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
    void CreateZeppelinCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
    void CreateTensorCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
    void CreateNODDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
    virtual void CreateDDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);

    CompartmentType m_CompartmentType;
    bool m_ModelWithFreeWaterComponent, m_ModelWithStationaryWaterComponent, m_ModelWithRestrictedWaterComponent;
    double m_FreeWaterProportionFixedValue, m_StationaryWaterProportionFixedValue, m_RestrictedWaterProportionFixedValue;
    unsigned int m_NumberOfCompartments;

    bool m_UseFixedWeights;
    bool m_UseConstrainedDiffusivity;
    bool m_UseConstrainedOrientationConcentration;
    bool m_UseConstrainedExtraAxonalFraction;
    bool m_UseConstrainedFreeWaterDiffusivity;
    bool m_UseConstrainedIRWDiffusivity;
    bool m_UseBoundedOptimization;

    bool m_UseCommonDiffusivities;
    bool m_UseCommonConcentrations;
    bool m_UseCommonExtraAxonalFractions;
    
    bool m_UserDefinedConcentrationBounds;
    double m_ConcentrationLowerBound;
    double m_ConcentrationUpperBound;

    double m_FreeWaterDiffusivity;
    double m_OrientationConcentration, m_ExtraAxonalFraction;
    double m_AxialDiffusivity;
    double m_RadialDiffusivity1, m_RadialDiffusivity2;
};

} // end namespace anima
