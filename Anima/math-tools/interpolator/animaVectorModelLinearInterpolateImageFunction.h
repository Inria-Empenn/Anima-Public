#pragma once

#include <itkInterpolateImageFunction.h>
#include <itkVariableLengthVector.h>

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
    typedef itk::VariableLengthVector <typename TInputImage::IOPixelType> VectorPixelType;

    /** Evaluate the function at a ContinuousIndex position
         *
         * Returns the linearly interpolated image intensity at a
         * specified point position. No bounds checking is done.
         * The point is assumed to lie within the image buffer.
         *
         * ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
    virtual OutputType EvaluateAtContinuousIndex(const ContinuousIndexType & index) const ITK_OVERRIDE;

protected:
    VectorModelLinearInterpolateImageFunction();
    virtual ~VectorModelLinearInterpolateImageFunction() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    virtual bool IsInsideBuffer(const ContinuousIndexType & index) const ITK_OVERRIDE
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

    template <class T> bool isZero(const itk::VariableLengthVector <T> &value) const
    {
        for (unsigned int i = 0;i < value.GetNumberOfElements();++i)
        {
            if (value[i] != 0)
                return false;
        }

        return true;
    }

    // Fake method for compilation purposes, should never go in there
    template <class T> bool isZero(T &data) const
    {
        itkExceptionMacro("Access to unauthorized method");
        return true;
    }

    //! Utility function to initialize output images pixel to zero for vector images
    template <class T> void InitializeZeroPixel(itk::VariableLengthVector <T> &zeroPixel) const
    {
        zeroPixel.Fill(0.0);
    }

    //! Utility function to initialize output images pixel to zero for all images except vector images
    template <class T> void InitializeZeroPixel(T &zeroPixel) const
    {
        zeroPixel = itk::NumericTraits <T>::ZeroValue();
    }

    //! Utility function to initialize output images pixel to zero for vector images
    template <class T1, class T2>
    void AddValue(itk::VariableLengthVector <T1> &input, double weight,
                  itk::VariableLengthVector <T2> &output) const
    {
        unsigned int vecDim = input.GetNumberOfElements();
        for (unsigned int i = 0;i < vecDim;++i)
            output[i] += weight * input[i];
    }

    //! Utility function to initialize output images pixel to zero for all images except vector images
    template <class T1, class T2> void AddValue(T1 &input, double weight, T2 &output) const
    {
        itkExceptionMacro("Acces to unauthorized method, check your implementation");
    }

private:
    VectorModelLinearInterpolateImageFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /** Number of neighbors used in the interpolation */
    static const unsigned long m_Neighbors;
};

} // end namespace itk

#include "animaVectorModelLinearInterpolateImageFunction.hxx"
