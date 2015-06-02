/**
 * SimuBlochSP-GRE
 * Version 0.3
 *
 * Description: A MRI simulator for the sequence of Spoiled Gradient Echo (possible called SPGR, FLASH, or T1-FFE on the scanner)
 *
 * Filename: SimuBlochSP-GRE.txx
 *
 * Libraries: This source file uses the open source library ITK 3.0
 * that is available through the following URI:
 * http://www.itk.org/
 * and the open source library TCLAP:
 * http://tclap.sourceforge.net/
 *
 *
 * @category   SimuBloch
 * @package    SimuBlochSP-GRE
 * @author     F. Cao <fang.cao@inria.fr>
 * @author     O. Commowick <Olivier.Commowick@inria.fr>
 * @copyright  INRIA Rennes - Bretagne Atlantique, VISAGES Team
 * @version    0.3
 *
 * Compile:    Install ITK
 *             Install TCLAP
 *             $Cmake
 *             $make
 *             (See readme for more information)
 *
 * Example:    $./SimuBlochSP-GRE  --t1 ./SampleData/T1.nii.gz --t2s ./SampleData/T2s.nii.gz --m0 ./SampleData/M0.nii.gz -o testing_SP-GRE_T1W.nii.gz --tr 35 --te 6 --fa 40
 *             (See readme for more information)
 *
 * Date:       Nov 16, 2012
 */

#pragma once
#include "animaSimuBlochSP-GRE.h"

#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template< class TImage>
SimuBlochSPGRE<TImage>::SimuBlochSPGRE()
{
    this->SetNumberOfRequiredInputs(3);
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputT1(const TImage* T1)
{
    SetInput(0, const_cast<TImage*>(T1));
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputT2s(const TImage* T2s)
{
    SetInput(1, const_cast<TImage*>(T2s));
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputM0(const TImage* M0)
{
    SetInput(2, const_cast<TImage*>(M0));
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputB1(const TImage* B1)
{
    this->SetInput(3, const_cast<TImage*>(B1));
}

template< class TImage>
void SimuBlochSPGRE< TImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{
    typename TImage::ConstPointer T1 = this->GetInput(0);
    typename TImage::ConstPointer T2s = this->GetInput(1);
    typename TImage::ConstPointer M0 = this->GetInput(2);
    typename TImage::ConstPointer B1;

    itk::ImageRegionIterator<TImage> outputIterator(this->GetOutput(), outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT1(T1, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2s(T2s, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorM0(M0, outputRegionForThread);

    itk::ImageRegionConstIterator<TImage> inputIteratorB1;
    bool b1DataPresent = (this->GetNumberOfIndexedInputs() == 4);
    if (b1DataPresent)
    {
        B1 = this->GetInput(3);
        inputIteratorB1 = itk::ImageRegionConstIterator<TImage> (B1, outputRegionForThread);
    }

    // For each voxel, calculate the signal of SP-GRE
    while(!outputIterator.IsAtEnd())
    {
        if ((inputIteratorT1.Get() <= 0) || ( inputIteratorT2s.Get() <= 0))
        {
            outputIterator.Set(0);

            ++inputIteratorT1;
            ++inputIteratorT2s;
            ++inputIteratorM0;
            if (b1DataPresent)
                ++inputIteratorB1;

            ++outputIterator;
            continue;
        }

        double b1Value = 1;
        if (b1DataPresent)
            b1Value = inputIteratorB1.Get();

        double t1Value = inputIteratorT1.Get();
        double t2sValue = inputIteratorT2s.Get();
        double m0Value = inputIteratorM0.Get();

        double sinFA = sin(itk::Math::pi * b1Value * m_FA / 180);
        double cosFA = cos(itk::Math::pi * b1Value * m_FA / 180);

        outputIterator.Set(m0Value * (1- exp(- m_TR / t1Value)) * sinFA / (1 - exp(- m_TR / t1Value) * cosFA ) * exp(- m_TE / t2sValue));

        ++inputIteratorT1;
        ++inputIteratorT2s;
        ++inputIteratorM0;
        if (b1DataPresent)
            ++inputIteratorB1;

        ++outputIterator;
    }

}
    
}// end of namespace anima
