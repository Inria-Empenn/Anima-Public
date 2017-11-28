#pragma once
#include "animaBaseTransformAgregator.h"

namespace anima
{

template <unsigned int NDimensions>
BaseTransformAgregator <NDimensions>::
BaseTransformAgregator()
{
    m_InputTransformType = TRANSLATION;
    m_OutputTransformType = TRANSLATION;
    m_Output = NULL;
    m_UpToDate = false;
    m_VerboseAgregation = true;

    m_InputTransforms.clear();
    m_InputOrigins.clear();
    m_Weights.clear();

    m_OrthogonalDirectionMatrix.Fill(0);
}

template <unsigned int NDimensions>
BaseTransformAgregator <NDimensions>::
~BaseTransformAgregator()
{
}

template <unsigned int NDimensions>
void
BaseTransformAgregator <NDimensions>::
SetInputTransformType(TRANSFORM_TYPE name)
{
    if (name == SVF)
    {
        std::cerr << "SVF input type is not yet supported, check your code" << std::endl;
        exit(-1);
    }

    if (name != m_InputTransformType)
        m_UpToDate = false;

    m_InputTransformType = name;
}

template <unsigned int NDimensions>
void
BaseTransformAgregator <NDimensions>::
SetOutputTransformType(TRANSFORM_TYPE name)
{
    if (name == DIRECTIONAL_AFFINE)
    {
        std::cerr << "Directional affine output type is not yet supported, check your code" << std::endl;
        exit(-1);
    }

    if (name != m_OutputTransformType)
        m_UpToDate = false;

    m_OutputTransformType = name;
}

template <unsigned int NDimensions>
typename BaseTransformAgregator <NDimensions>::BaseOutputTransformType *
BaseTransformAgregator <NDimensions>::
GetOutput()
{
    bool updateOk = true;
    if (!m_UpToDate)
        updateOk = this->Update();

    if (updateOk)
        return m_Output;
    else
        return NULL;
}

template <unsigned int NDimensions>
void
BaseTransformAgregator <NDimensions>::
SetOutput(BaseOutputTransformType *output)
{
    m_Output = output;
}

template <unsigned int NDimensions>
void
BaseTransformAgregator <NDimensions>::
SetInputTransforms(std::vector <BaseInputTransformPointer> &inputTransforms)
{
    m_InputTransforms = inputTransforms;

    m_UpToDate = false;
}


template <unsigned int NDimensions>
void
BaseTransformAgregator <NDimensions>::
SetCurrentLinearTransform(BaseInputTransformPointer &inputTransform)
{
    m_CurrentLinearTransform = inputTransform;

    m_UpToDate = false;
}

template <unsigned int NDimensions>
void
BaseTransformAgregator <NDimensions>::
SetOrthogonalDirectionMatrix(const MatrixType &inputTransform)
{
    m_OrthogonalDirectionMatrix = inputTransform;

    m_UpToDate = false;
}
} // end of namespace anima
