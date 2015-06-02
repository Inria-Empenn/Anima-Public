/**
 * SimuBlochGRE
 * Version 0.3
 *
 * Description: A MRI simulator for the sequence of Gradient Echo using Bloch equation
 *
 * Filename: SimuBlochGRE.txx
 *
 * Libraries: This source file uses the open source library ITK 3.0
 * that is available through the following URI:
 * http://www.itk.org/
 * and the open source library TCLAP:
 * http://tclap.sourceforge.net/
 *
 *
 * @category   SimuBloch
 * @package    SimuBlochGRE
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
 * Example:    $./SimuBlochGRE --t1 ./SampleData/T1.nii.gz --t2s ./SampleData/T2s.nii.gz --m0 ./SampleData/M0.nii.gz -o testing_GRE_T1W.nii.gz --tr 120 --te 8
 *             (See readme for more information)
 *
 * Date:       Nov 16, 2012
 */

#pragma once
#include "animaSimuBlochGRE.h"

#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{
    
template< class TImage>
SimuBlochGRE<TImage>::SimuBlochGRE()
{
    this->SetNumberOfRequiredInputs(3);
}

template< class TImage>
void SimuBlochGRE<TImage>::SetInputT1(const TImage* T1)
{
    SetInput(0, const_cast<TImage*>(T1));
}

template< class TImage>
void SimuBlochGRE<TImage>::SetInputT2s(const TImage* T2s)//changed for GRE
{
    SetInput(1, const_cast<TImage*>(T2s));
}

template< class TImage>
void SimuBlochGRE<TImage>::SetInputM0(const TImage* M0)
{
    SetInput(2, const_cast<TImage*>(M0));
}

template< class TImage>
void SimuBlochGRE< TImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{
    typename TImage::ConstPointer T1 = this->GetInput(0);
    typename TImage::ConstPointer T2s = this->GetInput(1); //changed for GRE
    typename TImage::ConstPointer M0 = this->GetInput(2);

    itk::ImageRegionIterator<TImage> outputIterator(this->GetOutput(), outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT1(T1, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2s(T2s, outputRegionForThread);  //changed for GRE
    itk::ImageRegionConstIterator<TImage> inputIteratorM0(M0, outputRegionForThread);

    // For each voxel, calculate the signal of Gradient Echo using bloch equation
    while(!outputIterator.IsAtEnd())
    {
        if ( (inputIteratorT1.Get() <= 0) || ( inputIteratorT2s.Get() <=0))  //changed for GRE
        {
            outputIterator.Set(0);
        }
        else
        {
            outputIterator.Set(inputIteratorM0.Get() * (1-exp(- m_TR / inputIteratorT1.Get())) * exp(- m_TE / inputIteratorT2s.Get()) ); //changed for GRE
        }

        ++inputIteratorT1;
        ++inputIteratorT2s; //changed for GRE
        ++inputIteratorM0;

        ++outputIterator;
    }

}
    
} // end of namespace anima
