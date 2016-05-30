#pragma once

#include <itkConfigure.h>
#include <itkInterpolateImageFunction.h>

namespace anima
{

template <class TInputImage, class TCoordRep = double>
class FasterLinearInterpolateImageFunction :
        public itk::InterpolateImageFunction<TInputImage,TCoordRep>
{
public:
    /** Standard class typedefs. */
    typedef FasterLinearInterpolateImageFunction                  Self;
    typedef itk::InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
    typedef itk::SmartPointer<Self>                              Pointer;
    typedef itk::SmartPointer<const Self>                        ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(FasterLinearInterpolateImageFunction, InterpolateImageFunction)

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** OutputType typedef support. */
    typedef typename Superclass::OutputType OutputType;

    /** InputImageType typedef support. */
    typedef typename Superclass::InputImageType InputImageType;

    /** InputPixelType typedef support. */
    typedef typename Superclass::InputPixelType InputPixelType;

    /** RealType typedef support. */
    typedef typename Superclass::RealType RealType;

    /** Dimension underlying input image. */
    itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

    /** Index typedef support. */
    typedef typename Superclass::IndexType      IndexType;
    typedef typename Superclass::IndexValueType IndexValueType;

    /** ContinuousIndex typedef support. */
    typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

    /** Evaluate the function at a ContinuousIndex position
         *
         * Returns the linearly interpolated image intensity at a
         * specified point position. No bounds checking is done.
         * The point is assume to lie within the image buffer.
         *
         * ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
    virtual OutputType EvaluateAtContinuousIndex(const ContinuousIndexType & index) const ITK_OVERRIDE;

    virtual bool IsInsideBuffer(const ContinuousIndexType & index) const ITK_OVERRIDE
    {
        for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
            /* Test for negative of a positive so we can catch NaN's. */
            if ( ! (index[j] >= this->m_StartIndex[j] &&
                    index[j] < this->m_EndIndex[j] ) )
            {
                return false;
            }
        }
        return true;
    }

protected:
    FasterLinearInterpolateImageFunction();
    virtual ~FasterLinearInterpolateImageFunction(){}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

private:
    FasterLinearInterpolateImageFunction( const Self& ); //purposely not implemented
    void operator=( const Self& ); //purposely not implemented

    /** Number of neighbors used in the interpolation */
    static const unsigned long  m_Neighbors;
};

} // end namespace of anima

# include "animaFasterLinearInterpolateImageFunction.hxx"
