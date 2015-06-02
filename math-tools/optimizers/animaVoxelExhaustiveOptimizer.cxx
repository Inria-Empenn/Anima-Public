
#include "animaVoxelExhaustiveOptimizer.h"
#include <itkCommand.h>
#include <itkEventObject.h>

namespace anima
{

/**
 * Constructor
 */
VoxelExhaustiveOptimizer
::VoxelExhaustiveOptimizer()
{
    m_CurrentIteration = 0;
    m_CurrentValue = 0;
    m_CurrentParameter = 0;
    m_CurrentIndex.Fill(0);
    m_Stop = false;
    m_NumberOfSteps.Fill(0);

    m_Geometry.set_size(3,3);
    m_Geometry.fill(0);
    m_Geometry.fill_diagonal(1);

    m_StopConditionDescription.str("");
}

const VoxelExhaustiveOptimizer::MeasureType& VoxelExhaustiveOptimizer
::GetCurrentCost() const
{
    if (this->m_Maximize)
        return this->GetMaximumMetricValue();
    else
        return this->GetMinimumMetricValue();
}

/**
 * Start walking
 */

void VoxelExhaustiveOptimizer::StartOptimization( void )
{
    this->StartWalking();

    if (this->m_Maximize)
    {
        this->m_CurrentValue = GetMaximumMetricValue();
        this->m_CurrentPosition = this->GetMaximumMetricValuePosition();
    }
    else
    {
        this->m_CurrentValue =  this->GetMinimumMetricValue();
        this->m_CurrentPosition = this->GetMinimumMetricValuePosition();
    }
}


void
VoxelExhaustiveOptimizer
::StartWalking( void )
{
    this->InvokeEvent( itk::StartEvent() );
    m_StopConditionDescription.str("");
    m_StopConditionDescription << this->GetNameOfClass() << ": Running";

    ParametersType initialPos = this->GetInitialPosition();
    m_MinimumMetricValuePosition = initialPos;
    m_MaximumMetricValuePosition = initialPos;

    MeasureType initialValue = this->GetValue( this->GetInitialPosition() );
    m_MaximumMetricValue = initialValue;
    m_MinimumMetricValue = initialValue;

    m_CurrentIteration          = 0;
    m_MaximumNumberOfIterations = 1;

    const unsigned int spaceDimension = this->GetInitialPosition().GetSize();

    for (unsigned int i=0; i< spaceDimension; i++)
    {
        m_MaximumNumberOfIterations *= (2 * m_NumberOfSteps[i] + 1);
    }

    m_CurrentIndex.SetSize(spaceDimension);
    m_CurrentIndex.Fill(0);

    ScalesType  scales = this->GetScales();

    // Make sure the scales have been set properly
    if (scales.size() != spaceDimension)
    {
        itkExceptionMacro(<< "The size of Scales is "
                          << scales.size()
                          << ", but the NumberOfParameters is "
                          << spaceDimension
                          << ".");
    }

    if (m_Geometry.rows() != spaceDimension)
    {
        itkExceptionMacro(<< "The size of the geometry matrix is "
                          << m_Geometry.rows()
                          << ", but the NumberOfParameters is "
                          << spaceDimension
                          << ".");
    }

    // Setup first grid position.
    ParametersType position( spaceDimension );
    for(unsigned int i=0; i<spaceDimension; i++)
        position[i] = this->GetInitialPosition()[i] - m_NumberOfSteps[i] * scales[i];

    position = m_Geometry * position;
    this->SetCurrentPosition( position );

    this->ResumeWalking();
}


/**
 * Resume the optimization
 */
void
VoxelExhaustiveOptimizer
::ResumeWalking( void )
{
    m_Stop = false;

    const unsigned int spaceDimension = this->GetInitialPosition().GetSize();
    ParametersType currentPosition(spaceDimension), newPosition(spaceDimension);

    while( !m_Stop )
    {
        currentPosition = this->GetCurrentPosition();

        if( m_Stop )
        {
            StopWalking();
            break;
        }

        m_CurrentValue = this->GetValue( currentPosition );

        if (m_CurrentValue > m_MaximumMetricValue)
        {
            m_MaximumMetricValue = m_CurrentValue;
            m_MaximumMetricValuePosition = currentPosition;
        }
        if (m_CurrentValue < m_MinimumMetricValue)
        {
            m_MinimumMetricValue = m_CurrentValue;
            m_MinimumMetricValuePosition = currentPosition;
        }

        if( m_Stop )
        {
            this->StopWalking();
            break;
        }

        m_StopConditionDescription.str("");
        m_StopConditionDescription << this->GetNameOfClass() << ": Running. ";
        m_StopConditionDescription << "@ index " << this->GetCurrentIndex() << " value is " << this->GetCurrentValue();

        this->InvokeEvent( itk::IterationEvent() );

        IncrementIndex(newPosition);
        newPosition = m_Geometry * newPosition;

        this->SetCurrentPosition(newPosition);

        m_CurrentIteration++;
    }
}


void
VoxelExhaustiveOptimizer
::StopWalking( void )
{
    m_Stop = true;
    this->InvokeEvent( itk::EndEvent() );
}


void
VoxelExhaustiveOptimizer
::IncrementIndex( ParametersType &newPosition )
{
    unsigned int idx = 0;
    const unsigned int  spaceDimension = m_CostFunction->GetNumberOfParameters();

    while( idx < spaceDimension )
    {
        m_CurrentIndex[idx]++;

        if( m_CurrentIndex[idx] > (2*m_NumberOfSteps[idx]))
        {
            m_CurrentIndex[idx]=0;
            idx++;
        }
        else
        {
            break;
        }
    }

    if( idx==spaceDimension )
    {
        m_Stop = true;
        m_StopConditionDescription.str("");
        m_StopConditionDescription << this->GetNameOfClass() << ": ";
        m_StopConditionDescription << "Completed sampling of parametric space of size " << spaceDimension;
    }

    for(unsigned int i=0; i<spaceDimension; i++)
        newPosition[i] = (m_CurrentIndex[i]-m_NumberOfSteps[i]) * this->GetScales()[i] + this->GetInitialPosition()[i];
}


const std::string
VoxelExhaustiveOptimizer
::GetStopConditionDescription() const
{
    return m_StopConditionDescription.str();
}

void
VoxelExhaustiveOptimizer
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "CurrentValue = " << m_CurrentValue << std::endl;
    os << indent << "NumberOfSteps = " << m_NumberOfSteps << std::endl;
    os << indent << "CurrentIteration = " << m_CurrentIteration << std::endl;
    os << indent << "Stop = " << m_Stop << std::endl;
    os << indent << "CurrentParameter = " << m_CurrentParameter << std::endl;
    os << indent << "CurrentIndex = " << m_CurrentIndex << std::endl;
    os << indent << "Geometry = " << m_Geometry << std::endl;
    os << indent << "MaximumNumberOfIterations = " << m_MaximumNumberOfIterations << std::endl;
    os << indent << "MaximumMetricValue = " << m_MaximumMetricValue << std::endl;
    os << indent << "MinimumMetricValue = " << m_MinimumMetricValue << std::endl;
    os << indent << "MinimumMetricValuePosition = " << m_MinimumMetricValuePosition << std::endl;
    os << indent << "MaximumMetricValuePosition = " << m_MaximumMetricValuePosition << std::endl;
}

} // end of namespace anima
