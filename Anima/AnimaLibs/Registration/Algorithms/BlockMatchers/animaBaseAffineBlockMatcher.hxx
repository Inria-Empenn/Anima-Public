#pragma once
#include <animaBaseAffineBlockMatcher.h>

/* Transforms */
#include <itkTranslationTransform.h>
#include <animaLogRigid3DTransform.h>
#include <animaSplitAffine3DTransform.h>
#include <animaDirectionScaleSkewTransform.h>
#include <animaNLOPTOptimizers.h>

namespace anima
{

template <typename TInputImageType>
BaseAffineBlockMatcher<TInputImageType>
::BaseAffineBlockMatcher()
{
    m_BlockTransformType = Translation;

    m_AngleMax = M_PI;
    m_TranslateMax = 10;
    m_ScaleMax = 3;

    m_AffineDirection = 1;
}

template <typename TInputImageType>
typename BaseAffineBlockMatcher<TInputImageType>::AgregatorType::TRANSFORM_TYPE
BaseAffineBlockMatcher<TInputImageType>
::GetAgregatorInputTransformType()
{
    switch (m_BlockTransformType)
    {
        case Translation:
            return AgregatorType::TRANSLATION;

        case Rigid:
            return AgregatorType::RIGID;

        case Affine:
        case Directional_Affine:
        default:
            return AgregatorType::AFFINE;
    }
}

template <typename TInputImageType>
typename BaseAffineBlockMatcher<TInputImageType>::BaseInputTransformPointer
BaseAffineBlockMatcher<TInputImageType>
::GetNewBlockTransform(PointType &blockCenter)
{
    BaseInputTransformPointer outputValue;

    switch(m_BlockTransformType)
    {
        case Translation:
        {
            typedef itk::TranslationTransform <double, InputImageType::ImageDimension> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();
            tr->SetIdentity();
            outputValue = tr;

            break;
        }

        case Rigid:
        {
            if (InputImageType::ImageDimension == 3)
            {
                typedef anima::LogRigid3DTransform<double> itkTransformType;
                typename itkTransformType::Pointer tr = itkTransformType::New();
                tr->SetIdentity();
                tr->SetCenter(blockCenter);

                outputValue = tr;
            }
            else
                std::cerr << "Only Rigid 3D transforms handled right now." << std::endl;
            break;
        }

        case Affine:
        {
            typedef anima::SplitAffine3DTransform <double> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();
            tr->SetIdentity();
            tr->SetCenter(blockCenter);

            outputValue = tr;
            break;
        }

        case Directional_Affine:
        default:
        {
            typedef anima::DirectionTransform <double> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();
            tr->SetIdentity();

            typename itkTransformType::HomogeneousMatrixType geometry;

            geometry.SetIdentity();
            for (unsigned int i = 0;i < 3;++i)
                for (unsigned int j = 0;j < 3;++j)
                    geometry(i,j) = this->GetReferenceImage()->GetDirection()(i,j) * this->GetReferenceImage()->GetSpacing()[j];

            for (unsigned int j = 0;j < 3;++j)
                geometry(j,TInputImageType::ImageDimension) = blockCenter[j];

            tr->SetUniqueDirection(m_AffineDirection);
            tr->SetGeometry(geometry, false);

            outputValue = tr;
            break;
        }
    }

    return outputValue;
}

template <typename TInputImageType>
void
BaseAffineBlockMatcher<TInputImageType>
::BlockMatchingSetup(MetricPointer &metric, unsigned int block)
{
    switch (m_BlockTransformType)
    {
        case Translation:
        {
            typedef itk::TranslationTransform <double, InputImageType::ImageDimension> itkTransformType;
            itkTransformType *tr = dynamic_cast <itkTransformType *> (this->GetBlockTransformPointer(block).GetPointer());
            tr->SetIdentity();
            break;
        }

        default:
        {
            typedef itk::MatrixOffsetTransformBase <double, InputImageType::ImageDimension> itkTransformType;
            itkTransformType *tr = dynamic_cast <itkTransformType *> (this->GetBlockTransformPointer(block).GetPointer());
            tr->SetIdentity();
            break;
        }
    }
}

template <typename TInputImageType>
void
BaseAffineBlockMatcher<TInputImageType>
::TransformDependantOptimizerSetup(OptimizerPointer &optimizer)
{
    if (this->GetOptimizerType() == Superclass::Exhaustive)
    {
        typedef anima::VoxelExhaustiveOptimizer LocalOptimizerType;
        LocalOptimizerType *tmpOpt = dynamic_cast <LocalOptimizerType *> (optimizer.GetPointer());

        LocalOptimizerType::ScalesType steps(InputImageType::ImageDimension);

        for (unsigned i = 0; i < steps.Size(); ++i)
            steps[i] = std::round(m_TranslateMax);

        tmpOpt->SetNumberOfSteps(steps);
        return;
    }

    typedef anima::NLOPTOptimizers LocalOptimizerType;
    LocalOptimizerType::ParametersType lowerBounds(this->GetBlockTransformPointer(0)->GetNumberOfParameters());
    LocalOptimizerType::ParametersType upperBounds(this->GetBlockTransformPointer(0)->GetNumberOfParameters());
    typename InputImageType::SpacingType fixedSpacing = this->GetReferenceImage()->GetSpacing();

    switch (m_BlockTransformType)
    {
        case Translation:
        {
            for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            {
                lowerBounds[i] = - m_TranslateMax * fixedSpacing[i];
                upperBounds[i] = m_TranslateMax * fixedSpacing[i];
            }

            break;
        }

        case Rigid:
        {
            for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            {
                lowerBounds[i] = - m_AngleMax * M_PI / 180.0;
                upperBounds[i] = m_AngleMax * M_PI / 180.0;

                lowerBounds[i+InputImageType::ImageDimension] = - m_TranslateMax * fixedSpacing[i];
                upperBounds[i+InputImageType::ImageDimension] = m_TranslateMax * fixedSpacing[i];
            }

            break;
        }

        case Affine:
        {
            // There are 12 parameters: 3 angles, 3 translations, 3 scales, 3 angles again
            for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            {
                // Angles
                lowerBounds[i] = - m_AngleMax * M_PI / 180.0;
                upperBounds[i] = m_AngleMax * M_PI / 180.0;

                // Translations
                lowerBounds[InputImageType::ImageDimension + i] = - m_TranslateMax * fixedSpacing[i];
                upperBounds[InputImageType::ImageDimension + i] = m_TranslateMax * fixedSpacing[i];

                // Scales
                lowerBounds[2 * InputImageType::ImageDimension + i] = - std::log (m_ScaleMax);
                upperBounds[2 * InputImageType::ImageDimension + i] = std::log (m_ScaleMax);

                // Second angles
                lowerBounds[3 * InputImageType::ImageDimension + i] = - m_AngleMax * M_PI / 180.0;
                upperBounds[3 * InputImageType::ImageDimension + i] = m_AngleMax * M_PI / 180.0;
            }

            break;
        }

        case Directional_Affine:
        default:
        {
            // There is 1 parameter: 1 translation in voxel coordinates
            lowerBounds[0] = - m_TranslateMax;
            upperBounds[0] = m_TranslateMax;

            break;
        }
    }

    LocalOptimizerType * tmpOpt = dynamic_cast <LocalOptimizerType *> (optimizer.GetPointer());
    tmpOpt->SetLowerBoundParameters(lowerBounds);
    tmpOpt->SetUpperBoundParameters(upperBounds);
}

} // end namespace anima
