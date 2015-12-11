#pragma once

#include <itkInterpolateImageFunction.h>

namespace anima
{

template <class TInputImage, class TCoordRep = double>
class VectorModelLinearInterpolateImageFunction :
        public itk::InterpolateImageFunction<TInputImage,TCoordRep>
{
public:
    /** Standard class typedefs. */
    typedef VectorModelLinearInterpolateImageFunction        Self;
    typedef itk::InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
    typedef itk::SmartPointer<Self>                              Pointer;
    typedef itk::SmartPointer<const Self>                        ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(VectorModelLinearInterpolateImageFunction,
                 InterpolateImageFunction);

    /** InputImageType typedef support. */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename TInputImage::PixelType     PixelType;
    typedef typename Superclass::RealType       RealType;

    /** Dimension underlying input image. */
    itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

    /** Index typedef support. */
    typedef typename Superclass::IndexType       IndexType;
    typedef typename Superclass::IndexValueType  IndexValueType;

    /** ContinuousIndex typedef support. */
    typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

    typedef typename Superclass::OutputType OutputType;

    /** Evaluate the function at a ContinuousIndex position
         *
         * Returns the linearly interpolated image intensity at a
         * specified point position. No bounds checking is done.
         * The point is assumed to lie within the image buffer.
         *
         * ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
    virtual OutputType EvaluateAtContinuousIndex(const ContinuousIndexType & index ) const;

protected:
    VectorModelLinearInterpolateImageFunction();
    virtual ~VectorModelLinearInterpolateImageFunction() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

    virtual bool IsInsideBuffer(const ContinuousIndexType & index) const
    {
        for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
            /* Test for negative of a positive so we can catch NaN's. */
            if ( ! (index[j] >= this->m_StartIndex[j] &&
                    index[j] <= this->m_EndIndex[j] ) )
            {
                return false;
            }
        }
        return true;
    }

    bool isZero(const PixelType &value) const
    {
        for (unsigned int i = 0;i < value.GetNumberOfElements();++i)
        {
            if (value[i] != 0)
            {
                return false;
            }
        }

        return true;
    }

private:
    VectorModelLinearInterpolateImageFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /** Number of neighbors used in the interpolation */
    static const unsigned long m_Neighbors;
};

} // end namespace itk

#include "animaVectorModelLinearInterpolateImageFunction.hxx"
