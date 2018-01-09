#pragma once

#include <itkSingleValuedNonLinearOptimizer.h>
#include <vnl/vnl_matrix.h>
#include "AnimaOptimizersExport.h"

namespace anima
{

class ANIMAOPTIMIZERS_EXPORT VoxelExhaustiveOptimizer :
public itk::SingleValuedNonLinearOptimizer
{
public:
    /** Standard "Self" typedef. */
    typedef VoxelExhaustiveOptimizer            Self;
    typedef itk::SingleValuedNonLinearOptimizer Superclass;
    typedef itk::SmartPointer<Self>             Pointer;
    typedef itk::SmartPointer<const Self>       ConstPointer;
    typedef vnl_matrix <double>            GeometryType;

    typedef itk::Array< unsigned long > StepsType;
    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(VoxelExhaustiveOptimizer, SingleValuedNonLinearOptimizer)

    virtual void StartOptimization() ITK_OVERRIDE;

    void StartWalking( void );
    void ResumeWalking( void );
    void StopWalking(void);

    itkSetMacro(NumberOfSteps, StepsType)

    // Geometry information
    itkSetMacro(Geometry, GeometryType)

    itkGetConstReferenceMacro(NumberOfSteps, StepsType)
    itkGetConstReferenceMacro(CurrentValue, MeasureType)
    itkGetConstReferenceMacro(MaximumMetricValue, MeasureType)
    itkGetConstReferenceMacro(MinimumMetricValue, MeasureType)
    itkGetConstReferenceMacro(MinimumMetricValuePosition, ParametersType)
    itkGetConstReferenceMacro(MaximumMetricValuePosition, ParametersType)
    itkGetConstReferenceMacro(CurrentIndex, ParametersType)
    itkGetConstReferenceMacro(MaximumNumberOfIterations, unsigned long)

    /** Get the reason for termination */
    const std::string GetStopConditionDescription() const ITK_OVERRIDE;

    /** Set if the Optimizer should Maximize the metric */
    itkSetMacro(Maximize, bool)
    itkBooleanMacro(Maximize)
    itkGetConstReferenceMacro(Maximize, bool)

    /** Return Current Value */

    const MeasureType& GetCurrentCost() const;

protected:
    VoxelExhaustiveOptimizer();
    virtual ~VoxelExhaustiveOptimizer() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    /** Advance to the next grid position. */
    void IncrementIndex(ParametersType &newPosition);


protected:
    MeasureType          m_CurrentValue;
    StepsType            m_NumberOfSteps;
    unsigned long        m_CurrentIteration;
    bool                 m_Stop;
    unsigned int         m_CurrentParameter;
    ParametersType       m_CurrentIndex;
    unsigned long        m_MaximumNumberOfIterations;
    MeasureType          m_MaximumMetricValue;
    MeasureType          m_MinimumMetricValue;
    ParametersType       m_MinimumMetricValuePosition;
    ParametersType       m_MaximumMetricValuePosition;

    GeometryType         m_Geometry;

private:
    VoxelExhaustiveOptimizer(const Self&); //purposely not implemented
    void operator=(const Self&);//purposely not implemented

    std::ostringstream m_StopConditionDescription;

    bool m_Maximize;
};

} // end of namespace anima
