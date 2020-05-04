#pragma once
#include "animaDWISimulatorFromDTIImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template <class PixelScalarType>
void
DWISimulatorFromDTIImageFilter<PixelScalarType>
::AddGradientDirection(unsigned int i, std::vector <double> &grad)
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
DWISimulatorFromDTIImageFilter<PixelScalarType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_GradientDirections.size());
}

template <class PixelScalarType>
void
DWISimulatorFromDTIImageFilter<PixelScalarType>
::BeforeThreadedGenerateData()
{
    if (this->GetInput()->GetNumberOfComponentsPerPixel() != m_NumberOfComponents)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Input image should have 6 dimensions...",ITK_LOCATION);

    Superclass::BeforeThreadedGenerateData();

    if (!m_S0Image)
    {
        m_S0Image = S0ImageType::New();
        m_S0Image->Initialize();
        m_S0Image->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        m_S0Image->SetSpacing (this->GetInput()->GetSpacing());
        m_S0Image->SetOrigin (this->GetInput()->GetOrigin());
        m_S0Image->SetDirection (this->GetInput()->GetDirection());
        m_S0Image->Allocate();

        m_S0Image->FillBuffer(m_S0Value);
    }
}

template <class PixelScalarType>
void
DWISimulatorFromDTIImageFilter<PixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;

    ImageIteratorType inIterator(this->GetInput(),outputRegionForThread);

    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    typedef itk::ImageRegionIterator <S0ImageType> S0ImageIteratorType;
    S0ImageIteratorType s0Iterator(m_S0Image,outputRegionForThread);

    typedef typename OutputImageType::PixelType OutputPixelType;

    unsigned int outputSize = m_GradientDirections.size();
    InputPixelType inVecTens(m_NumberOfComponents);
    OutputPixelType resVec(outputSize);

    while (!outIterator.IsAtEnd())
    {
        for (unsigned int i = 0;i < outputSize;++i)
            resVec[i] = 0;

        inVecTens = inIterator.Get();
        if (this->isZero(inVecTens))
        {
            outIterator.Set(resVec);

            ++inIterator;
            ++outIterator;
            ++s0Iterator;
            continue;
        }

        for (unsigned int i = 0;i < outputSize;++i)
        {
            unsigned int pos = 0;
            for (unsigned int j = 0;j < 3;++j)
                for (unsigned int k = 0;k <= j;++k)
                {
                    if (j != k)
                        resVec[i] += 2.0 * m_GradientDirections[i][j] * m_GradientDirections[i][k] * inVecTens[pos];
                    else
                        resVec[i] += m_GradientDirections[i][j] * m_GradientDirections[i][j] * inVecTens[pos];

                    ++pos;
                }

            resVec[i] *= - m_BValuesList[i];
            resVec[i] = s0Iterator.Get() * std::exp(resVec[i]);
        }

        outIterator.Set(resVec);

        ++inIterator;
        ++outIterator;
        ++s0Iterator;
    }
}

template <class PixelScalarType>
bool
DWISimulatorFromDTIImageFilter<PixelScalarType>
::isZero(InputPixelType &vec)
{
    unsigned int sizeVec = vec.GetSize();

    for (unsigned int i = 0;i < sizeVec;++i)
    {
        if (vec[i] != 0)
            return false;
    }

    return true;
}

template <class PixelScalarType>
typename DWISimulatorFromDTIImageFilter<PixelScalarType>::Image4DType *
DWISimulatorFromDTIImageFilter<PixelScalarType>
::GetOutputAs4DImage()
{
    this->Update();

    typename OutputImageType::Pointer outPtr = this->GetOutput();
    unsigned int ndim = outPtr->GetNumberOfComponentsPerPixel();

    typename Image4DType::RegionType region4d;
    typename InputImageType::RegionType region3d = outPtr->GetLargestPossibleRegion();
    typename Image4DType::SizeType size4d;
    typename Image4DType::SpacingType spacing4d;
    typename Image4DType::PointType origin4d;
    typename Image4DType::DirectionType direction4d;

    region4d.SetIndex(3,0);
    region4d.SetSize(3,ndim);
    size4d[3] = ndim;
    spacing4d[3] = 1;
    origin4d[3] = 0;
    direction4d(3,3) = 1;

    for (unsigned int i = 0;i < 3;++i)
    {
        region4d.SetIndex(i,region3d.GetIndex()[i]);
        region4d.SetSize(i,region3d.GetSize()[i]);

        size4d[i] = region3d.GetSize()[i];
        spacing4d[i] = outPtr->GetSpacing()[i];
        origin4d[i] = outPtr->GetOrigin()[i];
        for (unsigned int j = 0;j < 3;++j)
            direction4d(i,j) = (outPtr->GetDirection())(i,j);

        direction4d(i,3) = 0;
        direction4d(3,i) = 0;
    }

    m_Output4D = Image4DType::New();
    m_Output4D->Initialize();
    m_Output4D->SetRegions(region4d);
    m_Output4D->SetSpacing (spacing4d);
    m_Output4D->SetOrigin (origin4d);
    m_Output4D->SetDirection (direction4d);
    m_Output4D->Allocate();

    typedef itk::ImageRegionConstIterator <OutputImageType> OutImageIteratorType;
    typedef itk::ImageRegionIterator <Image4DType> Image4DIteratorType;

    for (unsigned int i = 0;i < ndim;++i)
    {
        region4d.SetIndex(3,i);
        region4d.SetSize(3,1);

        OutImageIteratorType outItr(outPtr,region3d);
        Image4DIteratorType fillItr(m_Output4D,region4d);

        while (!fillItr.IsAtEnd())
        {
            fillItr.Set(outItr.Get()[i]);
            ++fillItr;
            ++outItr;
        }
    }

    return m_Output4D;
}

} // end namespace anima
