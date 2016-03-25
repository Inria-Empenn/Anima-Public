/**
 * SimuBlochCoherentGRE
 * Version 0.3
 *
 * Description: A MRI simulator for the sequence of the Coherent or Partially Refocused (Rewound) Gradient Echo sequence (FISP, GRASS, FFE, FAST)
 *
 * Filename: SimuBlochCoherentGRE.txx
 *
 * Libraries: This source file uses the open source library ITK 3.0
 * that is available through the following URI:
 * http://www.itk.org/
 * and the open source library TCLAP:
 * http://tclap.sourceforge.net/
 *
 *
 * @category   SimuBloch
 * @package    SimuBlochCoherentGRE
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
 * Example:    $./SimuBlochCoherentGRE  --t1 ./SampleData/T1.nii.gz --t2 ./SampleData/T2.nii.gz --t2s ./SampleData/T2s.nii.gz --m0 ./SampleData/M0.nii.gz -o testing_CoherentGRE_T2sW.nii.gz --tr 40 --te 15 --fa 25
 *             (See readme for more information)
 *
 * Date:       Nov 16, 2012
 */

#pragma once
#include "animaSimuBlochCoherentGRE.h"

#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template< class TImage>
SimuBlochCoherentGRE<TImage>::SimuBlochCoherentGRE()
{
    this->SetNumberOfRequiredInputs(3);
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputT1(const TImage* T1)
{
    SetInput(0, const_cast<TImage*>(T1));
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputT2s(const TImage* T2s)
{
    SetInput(1, const_cast<TImage*>(T2s));
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputM0(const TImage* M0)
{
    SetInput(2, const_cast<TImage*>(M0));
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputT2(const TImage* T2)//changed for CoherentGRE
{
    SetInput(3, const_cast<TImage*>(T2));//changed for CoherentGRE
}

template< class TImage>//changed for CoherentGRE

void SimuBlochCoherentGRE< TImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{
    typename TImage::ConstPointer T1 = this->GetInput(0);
    typename TImage::ConstPointer T2s = this->GetInput(1);
    typename TImage::ConstPointer M0 = this->GetInput(2);
    typename TImage::ConstPointer T2 = this->GetInput(3);//changed for CoherentGRE

    itk::ImageRegionIterator<TImage> outputIterator(this->GetOutput(), outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT1(T1, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2s(T2s, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorM0(M0, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2(T2, outputRegionForThread);//changed for CoherentGRE

    float sinFA= sin(M_PI * m_FA / 180);
    float cosFA = cos(M_PI * m_FA / 180);

    // For each voxel, calculate the signal of CoherentGRE
    while(!outputIterator.IsAtEnd())
    {
        if ( (inputIteratorT1.Get() <= 0) || ( inputIteratorT2s.Get() <=0) || ( inputIteratorT2.Get() <=0))//changed for CoherentGRE
        {
            outputIterator.Set(0);
        }
        else
        {
            //      outputIterator.Set(inputIteratorM0.Get() * (1- exp(- m_TR / inputIteratorT1.Get())) * sinFA / (1 - exp(- m_TR / inputIteratorT1.Get()) * cosFA  - exp( - m_TR /  inputIteratorT2.Get() ) * (exp( - m_TR / inputIteratorT1.Get() ) - cosFA) ) * exp(- m_TE / inputIteratorT2s.Get()) ); //changed for CoherentGRE//changed for CoherentGRE
            outputIterator.Set(inputIteratorM0.Get() * sinFA / (1 + inputIteratorT1.Get() / inputIteratorT2.Get() - cosFA * ( inputIteratorT1.Get()/ inputIteratorT2.Get() -1 ) ) * exp(- m_TE / inputIteratorT2s.Get()) ); //changed for CoherentGRE//changed for CoherentGRE
        }

        ++inputIteratorT1;
        ++inputIteratorT2s;
        ++inputIteratorM0;
        ++inputIteratorT2;//changed for CoherentGRE

        ++outputIterator;
    }
}
    
}// end of namespace anima
