#pragma once
#include "animaDTIEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include <nlopt.hpp>

#include <animaVectorOperations.h>
#include <animaBaseTensorTools.h>

#include <boost/lexical_cast.hpp>

namespace anima
{

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::AddGradientDirection(unsigned int i, vnl_vector_fixed<double,3> &grad)
{
    if (i == m_GradientDirections.size())
        m_GradientDirections.push_back(grad);
    else if (i > m_GradientDirections.size())
        std::cerr << "Trying to add a direction not contiguous... Add directions contiguously (0,1,2,3,...)..." << std::endl;
    else
        m_GradientDirections[i] = grad;
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;
    
    typedef itk::ImageRegionConstIterator <InputImageType> B0IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    unsigned int firstB0Index = 0;
    while (m_BValuesList[firstB0Index] > 10)
        ++firstB0Index;
    
    B0IteratorType b0Itr(this->GetInput(firstB0Index),this->GetOutput()->GetLargestPossibleRegion());

    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    maskImage->SetSpacing (this->GetInput(0)->GetSpacing());
    maskImage->SetOrigin (this->GetInput(0)->GetOrigin());
    maskImage->SetDirection (this->GetInput(0)->GetDirection());
    maskImage->Allocate();

    MaskIteratorType maskItr (maskImage,this->GetOutput()->GetLargestPossibleRegion());
    while (!b0Itr.IsAtEnd())
    {
        if (b0Itr.Get() > m_B0Threshold)
            maskItr.Set(1);
        else
            maskItr.Set(0);

        ++b0Itr;
        ++maskItr;
    }

    this->SetComputationMask(maskImage);
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfComponents);
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::BeforeThreadedGenerateData()
{
    if (m_BValuesList.size() != this->GetNumberOfIndexedInputs())
    {
        std::string error("There should be the same number of input images and input b-values... ");
        error += boost::lexical_cast <std::string> (m_BValuesList.size());
        error += " ";
        error += boost::lexical_cast <std::string> (this->GetNumberOfIndexedInputs());

        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    itk::ImageRegionIterator <OutputImageType> fillOut(this->GetOutput(),this->GetOutput()->GetLargestPossibleRegion());
    typename OutputImageType::PixelType emptyModelVec(m_NumberOfComponents);
    emptyModelVec.Fill(0);

    while (!fillOut.IsAtEnd())
    {
        fillOut.Set(emptyModelVec);
        ++fillOut;
    }

    m_EstimatedB0Image = OutputB0ImageType::New();
    m_EstimatedB0Image->Initialize();
    m_EstimatedB0Image->SetOrigin(this->GetOutput()->GetOrigin());
    m_EstimatedB0Image->SetSpacing(this->GetOutput()->GetSpacing());
    m_EstimatedB0Image->SetDirection(this->GetOutput()->GetDirection());
    m_EstimatedB0Image->SetRegions(this->GetOutput()->GetLargestPossibleRegion());

    m_EstimatedB0Image->Allocate();
    m_EstimatedB0Image->FillBuffer(0.0);

    Superclass::BeforeThreadedGenerateData();
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionConstIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskIterator(this->GetComputationMask(),outputRegionForThread);
    
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    typedef itk::ImageRegionIteratorWithIndex <OutputB0ImageType> OutB0ImageIteratorType;
    OutB0ImageIteratorType outB0Iterator(m_EstimatedB0Image,outputRegionForThread);

    typedef typename OutputImageType::PixelType OutputPixelType;
    std::vector <double> dwi (numInputs,0);
    OutputPixelType resVec(m_NumberOfComponents);

    vnl_matrix <double> tmpTensor(3,3);

    std::vector <double> predictedValues(numInputs,0);
    OptimizationDataStructure data;

    while (!outIterator.IsAtEnd())
    {
        resVec.Fill(0.0);
        
        if (maskIterator.Get() == 0)
        {
            outIterator.Set(resVec);
            outB0Iterator.Set(0);
            
            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];
            
            ++maskIterator;
            ++outIterator;
            ++outB0Iterator;
            
            continue;
        }

        for (unsigned int i = 0;i < numInputs;++i)
            dwi[i] = inIterators[i].Get();

        // NLOPT optimization
        nlopt::opt opt(nlopt::LN_BOBYQA, m_NumberOfComponents);

        std::vector <double> lowerBounds(m_NumberOfComponents, 0);
        for (unsigned int i = 0;i < 3;++i)
            lowerBounds[i + 3] = -1.0e-3;

        opt.set_lower_bounds(lowerBounds);

        std::vector <double> upperBounds(m_NumberOfComponents, 5.0e-3);
        for (unsigned int i = 0;i < 3;++i)
            upperBounds[i + 3] = 1.0e-3;

        opt.set_upper_bounds(upperBounds);
        opt.set_xtol_rel(1e-4);

        std::vector <double> optimizedValue(m_NumberOfComponents, 1.0e-3);
        for (unsigned int i = 0;i < 3;++i)
            upperBounds[i + 3] = 0.0;

        double minf;

        data.filter = this;
        data.dwi = dwi;
        data.predictedValues = predictedValues;

        opt.set_min_objective(OptimizationFunction, &data);

        try
        {
            opt.optimize(optimizedValue, minf);
        }
        catch(nlopt::roundoff_limited& e)
        {
            itkExceptionMacro("NLOPT optimization error");
        }

        // Tensor is column first
        resVec[0] = optimizedValue[0];
        resVec[2] = optimizedValue[1];
        resVec[5] = optimizedValue[2];
        resVec[1] = optimizedValue[3];
        resVec[3] = optimizedValue[4];
        resVec[4] = optimizedValue[5];

        double outB0Value = 0;
        double normConstant = 0;

        for (unsigned int i = 0;i < numInputs;++i)
        {
            double predictedValue = 0;
            double bValue = m_BValuesList[i];

            if (bValue == 0)
                predictedValue = 1;
            else
            {
                for (unsigned int j = 0;j < 3;++j)
                    predictedValue += optimizedValue[j] * m_GradientDirections[i][j] * m_GradientDirections[i][j];

                predictedValue += 2 * optimizedValue[3] * m_GradientDirections[i][0] * m_GradientDirections[i][1];
                predictedValue += 2 * optimizedValue[4] * m_GradientDirections[i][0] * m_GradientDirections[i][2];
                predictedValue += 2 * optimizedValue[5] * m_GradientDirections[i][1] * m_GradientDirections[i][2];

                predictedValue = std::exp(- bValue * predictedValue);
            }

            outB0Value += predictedValue * dwi[i];
            normConstant += predictedValue * predictedValue;
        }

        outB0Value /= normConstant;

        anima::GetTensorFromVectorRepresentation(resVec,tmpTensor,3);
        vnl_symmetric_eigensystem <double> tmpEigs(tmpTensor);

        bool isTensorOk = true;

        if (m_RemoveDegeneratedTensors)
        {
            for (unsigned int i = 0;i < 3;++i)
            {
                if (tmpEigs.D[i] <= 0)
                {
                    isTensorOk = false;
                    break;
                }
            }
        }

        if (isTensorOk)
        {
            outIterator.Set(resVec);
            outB0Iterator.Set(outB0Value);
        }

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        ++maskIterator;
        ++outIterator;
        ++outB0Iterator;
    }
}

template <class PixelScalarType>
double
DTIEstimationImageFilter<PixelScalarType>
::OptimizationFunction(const std::vector<double> &x, std::vector<double> &grad, void *func_data)
{
    OptimizationDataStructure *data = static_cast <OptimizationDataStructure *> (func_data);

    return data->filter->ComputeCostAtPosition(x,data->dwi,data->predictedValues);
}

template <class PixelScalarType>
double
DTIEstimationImageFilter<PixelScalarType>
::ComputeCostAtPosition(const std::vector<double> &x, const std::vector <double> &observedData,
                        std::vector <double> &predictedValues)
{
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        predictedValues[i] = 0;
        double bValue = m_BValuesList[i];
        if (bValue == 0)
        {
            predictedValues[i] = 1;
            continue;
        }

        for (unsigned int j = 0;j < 3;++j)
            predictedValues[i] += x[j] * m_GradientDirections[i][j] * m_GradientDirections[i][j];

        predictedValues[i] += 2 * x[3] * m_GradientDirections[i][0] * m_GradientDirections[i][1];
        predictedValues[i] += 2 * x[4] * m_GradientDirections[i][0] * m_GradientDirections[i][2];
        predictedValues[i] += 2 * x[5] * m_GradientDirections[i][1] * m_GradientDirections[i][2];

        predictedValues[i] = std::exp(- bValue * predictedValues[i]);
    }

    double b0Val = 0;
    double normConstant = 0;

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        b0Val += predictedValues[i] * observedData[i];
        normConstant += predictedValues[i] * predictedValues[i];
    }

    b0Val /= normConstant;

    double costValue = 0;
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        double predVal = b0Val * predictedValues[i];
        costValue += (predVal - observedData[i]) * (predVal - observedData[i]);
    }

    return costValue;
}

} // end of namespace anima
