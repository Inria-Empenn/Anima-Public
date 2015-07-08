#pragma once

#include <itkInPlaceImageFilter.h>
#include <itkNumericTraits.h>
#include <itkImageRegionSplitterDirection.h>
#include <itkVector.h>
#include <itkVariableLengthVector.h>

namespace anima
{
template <typename TInputImage, typename TOutputImage=TInputImage>
class RecursiveLineYvvGaussianImageFilter :
public itk::InPlaceImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef RecursiveLineYvvGaussianImageFilter                  Self;
    typedef itk::InPlaceImageFilter<TInputImage,TOutputImage>   Superclass;
    typedef itk::SmartPointer<Self>                             Pointer;
    typedef itk::SmartPointer<const Self>                       ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Type macro that defines a name for this class. */
    itkTypeMacro(RecursiveLineYvvGaussianImageFilter, InPlaceImageFilter);

    /** Smart pointer typedef support.  */
    typedef typename TInputImage::Pointer       InputImagePointer;
    typedef typename TInputImage::ConstPointer  InputImageConstPointer;

    /** Real type to be used in internal computations. RealType in general is
     * templated over the pixel type. (For example for vector or tensor pixels,
     * RealType is a vector or a tensor of doubles.) ScalarRealType is a type
     * meant for scalars.
     */
    typedef typename TInputImage::PixelType                        InputPixelType;
    typedef typename itk::NumericTraits<InputPixelType>::RealType       RealType;
    typedef typename itk::NumericTraits<InputPixelType>::ScalarRealType ScalarRealType;

    typedef typename TOutputImage::RegionType                      OutputImageRegionType;

    /** Type of the input image */
    typedef TInputImage      InputImageType;

    /** Type of the output image */
    typedef TOutputImage      OutputImageType;

    /** Get the direction in which the filter is to be applied. */
    itkGetConstMacro(Direction, unsigned int);

    /** Set the direction in which the filter is to be applied. */
    itkSetMacro(Direction, unsigned int);

    /** Set Input Image. */
    void SetInputImage( const TInputImage * );

    /** Get Input Image. */
    const TInputImage * GetInputImage( void );

    /** Set/Get the flag for normalizing the gaussian over scale space.
     When this flag is ON the filter will be normalized in such a way
     that larger sigmas will not result in the image fading away.

     \f[
     \frac{ 1 }{ \sqrt{ 2 \pi } };
     \f]

     When the flag is OFF the normalization will conserve contant the
     integral of the image intensity.
     \f[
     \frac{ 1 }{ \sigma  \sqrt{ 2 \pi } };
     \f]
     For analyzing an image across Scale Space you want to enable
     this flag.  It is disabled by default.  */
    itkSetMacro( NormalizeAcrossScale, bool );
    itkGetConstMacro( NormalizeAcrossScale, bool );

    /** Set/Get the Sigma, measured in world coordinates, of the Gaussian
     * kernel.  The default is 1.0.  */
    itkGetConstMacro( Sigma, ScalarRealType );
    itkSetMacro( Sigma, ScalarRealType );

protected:
    RecursiveLineYvvGaussianImageFilter();
    virtual ~RecursiveLineYvvGaussianImageFilter() {};
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

    /** GenerateData (apply) the filter. */
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId );

    virtual const itk::ImageRegionSplitterBase* GetImageRegionSplitter(void) const;

    /** RecursiveLineYvvGaussianImageFilter needs all of the input only in the
     *  "Direction" dimension. Therefore we enlarge the output's
     *  RequestedRegion to this. Then the superclass's
     *  GenerateInputRequestedRegion method will copy the output region
     *  to the input.
     *
     * \sa ImageToImageFilter::GenerateInputRequestedRegion()
     */
    void EnlargeOutputRequestedRegion(itk::DataObject *output);

    /** Set up the coefficients of the filter to approximate a specific kernel.
     * Typically it can be used to approximate a Gaussian or one of its
     * derivatives. Parameter is the spacing along the dimension to
     * filter. */
    virtual void SetUp(ScalarRealType spacing);

    /** Apply the Recursive Filter to an array of data.  This method is called
     * for each line of the volume. Parameter "scratch" is a scratch
     * area used for internal computations that is the same size as the
     * parameters "outs" and "data". The scratch area must be allocated
     * outside of this routine (this avoids memory allocation and
     * deallocation in the inner loop of the overall algorithm. */
    void FilterDataArray(RealType *outs, const RealType *data, unsigned int ln,
                         RealType &sV0, RealType &sV1, RealType &sV2);

protected:
    /** Utility function to compute causal part for any type */
    template <class T>
    inline void ComputeCausalPart(T &out, const T &data, T &V0, T &V1, T &V2)
    {
        out = data + V0 * m_B1 + V1 * m_B2 + V2 * m_B3;
        V2 = V1;V1 = V0;V0 = out;
    }

    /** Utility function to compute causal part for ITK vector type */
    template <class T>
    inline void ComputeCausalPart(itk::VariableLengthVector<T> &out, const itk::VariableLengthVector<T> &data,
                                  itk::VariableLengthVector<T> &V0, itk::VariableLengthVector<T> &V1,
                                  itk::VariableLengthVector<T> &V2)
    {
        unsigned int vSize = data.GetSize();
        if (out.GetSize() != vSize)
            out.SetSize(vSize);

        for (unsigned int i = 0;i < vSize;++i)
        {
            out[i] = data[i] + V0[i] * m_B1 + V1[i] * m_B2 + V2[i] * m_B3;
            V2[i] = V1[i];V1[i] = V0[i];V0[i] = out[i];
        }
    }

    /** Utility function to compute causal part for any type */
    template <class T>
    inline void ComputeAntiCausalPart(T &out, T &data, T &V0, T &V1, T &V2)
    {
        out = data * m_B + V0 * m_B1 + V1 * m_B2 + V2 * m_B3;
        V2 = V1;V1 = V0;V0 = out;
    }

    /** Utility function to compute causal part for ITK vector type */
    template <class T>
    inline void ComputeAntiCausalPart(itk::VariableLengthVector<T> &out, itk::VariableLengthVector<T> &data,
                                      itk::VariableLengthVector<T> &V0, itk::VariableLengthVector<T> &V1,
                                      itk::VariableLengthVector<T> &V2)
    {
        unsigned int vSize = data.GetSize();
        for (unsigned int i = 0;i < vSize;++i)
        {
            out[i] = data[i] * m_B + V0[i] * m_B1 + V1[i] * m_B2 + V2[i] * m_B3;
            V2[i] = V1[i];V1[i] = V0[i];V0[i] = out[i];
        }
    }

    template <class T>
    inline void ComputeCausalBase(const T &data, T &V0, T &V1, T &V2)
    {
        V0 = data / (1.0 - m_B1 - m_B2 - m_B3);V1 = V0;V2 = V1;
    }

    template <class T>
    inline void ComputeCausalBase(const itk::VariableLengthVector<T> &data, itk::VariableLengthVector<T> &V0,
                                  itk::VariableLengthVector<T> &V1, itk::VariableLengthVector<T> &V2)
    {
        unsigned int vSize = data.GetSize();
        if (V0.GetSize() != vSize)
            V0.SetSize(vSize);
        if (V1.GetSize() != vSize)
            V1.SetSize(vSize);
        if (V2.GetSize() != vSize)
            V2.SetSize(vSize);

        for (unsigned int i = 0;i < vSize;++i)
        {
            V0[i] = data[i] / (1.0 - m_B1 - m_B2 - m_B3);V1[i] = V0[i];V2[i] = V1[i];
        }
    }

    template <class T>
    inline void ComputeAntiCausalBase(const T &data, T *outs, T &V0, T &V1, T &V2, unsigned int ln)
    {
        // Handle outside values according to Triggs and Sdika
        const T u_p = data / (1.0 - m_B1 - m_B2 - m_B3);
        const T v_p = u_p / (1.0 - m_B1 - m_B2 - m_B3);

        V0 = v_p;
        V1 = v_p;
        V2 = v_p;

        for (unsigned int i = 0;i < 3;++i)
        {
            V0 += (outs[ln - 1 - i] - u_p) * m_MMatrix(0,i);
            V1 += (outs[ln - 1 - i] - u_p) * m_MMatrix(1,i);
            V2 += (outs[ln - 1 - i] - u_p) * m_MMatrix(2,i);
        }

        // This was not in the 2006 Triggs paper but sounds quite logical since m_B is not one
        V0 *= m_B;
        V1 *= m_B;
        V2 *= m_B;
    }

    template <class T>
    inline void ComputeAntiCausalBase(const itk::VariableLengthVector<T> &data, itk::VariableLengthVector<T> *outs,
                                      itk::VariableLengthVector<T> &V0, itk::VariableLengthVector<T> &V1,
                                      itk::VariableLengthVector<T> &V2, unsigned int ln)
    {
        unsigned int vSize = data.GetSize();
        // Handle outside values according to Triggs and Sdika
        double factor = 1.0 / (1.0 - m_B1 - m_B2 - m_B3);

        for (unsigned int i = 0;i < vSize;++i)
        {
            V0[i] = data[i] * factor * factor;
            V1[i] = data[i] * factor * factor;
            V2[i] = data[i] * factor * factor;
        }

        for (unsigned int i = 0;i < 3;++i)
            for (unsigned int j = 0;j < vSize;++j)
        {
            V0[j] += (outs[ln - 1 - i][j] - data[j] * factor) * m_MMatrix(0,i);
            V1[j] += (outs[ln - 1 - i][j] - data[j] * factor) * m_MMatrix(1,i);
            V2[j] += (outs[ln - 1 - i][j] - data[j] * factor) * m_MMatrix(2,i);
        }

        // This was not in the 2006 Triggs paper but sounds quite logical since m_B is not one
        for (unsigned int i = 0;i < vSize;++i)
        {
            V0[i] *= m_B;
            V1[i] *= m_B;
            V2[i] *= m_B;
        }
    }

    /** Causal and anti-causal coefficients that multiply the input data. These are already divided by B0 */
    ScalarRealType m_B1;
    ScalarRealType m_B2;
    ScalarRealType m_B3;
    ScalarRealType m_B;

    // Initialization matrix for anti-causal pass
    vnl_matrix <ScalarRealType> m_MMatrix;

private:
    RecursiveLineYvvGaussianImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /** Direction in which the filter is to be applied
     * this should be in the range [0,ImageDimension-1]. */
    unsigned int m_Direction;

    /** Sigma of the gaussian kernel. */
    ScalarRealType m_Sigma;

    /** Normalize the image across scale space */
    bool m_NormalizeAcrossScale;

    itk::ImageRegionSplitterDirection::Pointer m_ImageRegionSplitter;
};


} // end of namespace anima

#include "animaRecursiveLineYvvGaussianImageFilter.hxx"
