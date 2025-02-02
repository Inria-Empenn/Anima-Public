#pragma once

#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>
#include <itkImageRegionSplitterDirection.h>
#include <itkVector.h>

namespace anima
{

template <typename TInputImage, typename TDiffusionScalarImage=TInputImage, typename TOutputImage=TInputImage>
class InhomogeneousAOSLineDiffusionImageFilter :
        public itk::ImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef InhomogeneousAOSLineDiffusionImageFilter       Self;
    typedef itk::ImageToImageFilter<TInputImage,TOutputImage>   Superclass;
    typedef itk::SmartPointer<Self>                             Pointer;
    typedef itk::SmartPointer<const Self>                       ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Type macro that defines a name for this class. */
    itkTypeMacro(InhomogeneousAOSLineDiffusionImageFilter, ImageToImageFilter)

    /** Smart pointer typedef support.  */
    typedef typename TInputImage::Pointer       InputImagePointer;
    typedef typename TInputImage::ConstPointer  InputImageConstPointer;

    /** Real type to be used in internal computations. RealType in general is
         * templated over the pixel type. (For example for vector or tensor pixels,
         * RealType is a vector or a tensor of doubles.) ScalarRealType is a type
         * meant for scalars.
         */
    typedef typename TInputImage::PixelType InputPixelType;
    typedef typename itk::NumericTraits<InputPixelType>::RealType RealType;
    typedef typename itk::NumericTraits<InputPixelType>::ScalarRealType ScalarRealType;

    typedef typename TOutputImage::RegionType OutputImageRegionType;

    /** Type of the input image */
    typedef TInputImage InputImageType;

    typedef TDiffusionScalarImage DiffusionScalarsImageType;
    typedef typename DiffusionScalarsImageType::Pointer DiffusionScalarsImagePointer;
    typedef TOutputImage OutputImageType;

    /** Get the direction in which the filter is to be applied. */
    itkGetConstMacro(Direction, unsigned int)

    /** Set the direction in which the filter is to be applied. */
    itkSetMacro(Direction, unsigned int)

    /** Set Input Image. */
    void SetInputImage( const TInputImage * );

    /** Get Input Image. */
    const TInputImage * GetInputImage();

    itkSetMacro(TimeStep, double)
    void SetDiffusionScalarsImage (DiffusionScalarsImageType *data) {m_DiffusionScalarsImage = data;}

protected:
    InhomogeneousAOSLineDiffusionImageFilter();
    virtual ~InhomogeneousAOSLineDiffusionImageFilter() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    /** GenerateData (apply) the filter. */
    void GenerateData() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) ITK_OVERRIDE;

    /** InhomogeneousAOSLineDiffusionImageFilter needs all of the input only in the
         *  "Direction" dimension. Therefore we enlarge the output's
         *  RequestedRegion to this. Then the superclass's
         *  GenerateInputRequestedRegion method will copy the output region
         *  to the input.
         *
         * \sa ImageToImageFilter::GenerateInputRequestedRegion()
         */
    void EnlargeOutputRequestedRegion(itk::DataObject *output) ITK_OVERRIDE;

    /** Apply the filter to an array of data.  This method is called
         * for each line of the volume. Parameter "scratch" is a scratch
         * area used for internal computations that is the same size as the
         * parameters "outs" and "data". l_coefs, m_coefs and r_coefs are temporary variables
         * to hold LR decomposition results. */
    void FilterDataArray(std::vector <RealType> &outs, const std::vector <RealType> &data,
                         std::vector <RealType> &diffs, std::vector <RealType> &scratch,
                         std::vector <RealType> &l_coefs, std::vector <RealType> &m_coefs,
                         std::vector <RealType> &r_coefs);

    template <typename T1, typename T2>
    void
    InitializeLeftSideCoefs(T1 &outL0, T1 &outM0, T1 &outR0,
                            const T2 &delta, const T1 &diffs0,
                            const T1 &diffs1)
    {
        outL0 = 0;
        outM0 = 1.0 + 2 * delta * diffs0 / m_SquareGradientDelta;
        outR0 = - delta * (5 * diffs0 - diffs1) / (4 * m_SquareGradientDelta);
    }

    template <typename T1, typename T2, unsigned int TDimension>
    void
    InitializeLeftSideCoefs(itk::Vector <T1,TDimension> &outL0,
                            itk::Vector <T1,TDimension> &outM0,
                            itk::Vector <T1,TDimension> &outR0,
                            const T2 &delta, const itk::Vector <T1,TDimension> &diffs0,
                            const itk::Vector <T1,TDimension> &diffs1)
    {
        for (unsigned int i = 0;i < TDimension;++i)
        {
            outL0[i] = 0;
            outM0[i] = 1.0 + 2 * delta * diffs0[i] / m_SquareGradientDelta;
            outR0[i] = - delta * (5 * diffs0[i] - diffs1[i]) / (4 * m_SquareGradientDelta);
        }
    }

    template <typename T1, typename T2>
    void
    InitializeIthCoefs(T1 &outLi, T1 &outMi, T1 &outRi, const T2 &delta,
                       const T1 &inPrevM, const T1 &inPrevR,
                       const T1 &diffsPrev, const T1 &diffs,
                       const T1 &diffsNext)
    {
        outLi = 1.0 / inPrevM;
        outLi *= - delta * (diffs + (diffsNext - diffsPrev) / 4.0) / m_SquareGradientDelta;

        outMi = 1.0 + 2 * delta * diffs / m_SquareGradientDelta;
        outMi -= outLi * inPrevR;

        outRi = - delta * (diffs - (diffsNext - diffsPrev) / 4.0) / m_SquareGradientDelta;
    }

    template <typename T1, typename T2, unsigned int TDimension>
    void
    InitializeIthCoefs(itk::Vector <T1,TDimension> &outLi,
                       itk::Vector <T1,TDimension> &outMi,
                       itk::Vector <T1,TDimension> &outRi, const T2 &delta,
                       const itk::Vector <T1,TDimension> &inPrevM,
                       const itk::Vector <T1,TDimension> &inPrevR,
                       const itk::Vector <T1,TDimension> &diffsPrev,
                       const itk::Vector <T1,TDimension> &diffs,
                       const itk::Vector <T1,TDimension> &diffsNext)
    {
        for (unsigned int i = 0;i < TDimension;++i)
        {
            outLi[i] = 1.0 / inPrevM[i];
            outLi[i] *= - delta * (diffs[i] + (diffsNext[i] - diffsPrev[i]) / 4.0) / m_SquareGradientDelta;

            outMi[i] = 1.0 + 2 * delta * diffs[i] / m_SquareGradientDelta;
            outMi[i] -= outLi[i] * inPrevR[i];

            outRi[i] = - delta * (diffs[i] - (diffsNext[i] - diffsPrev[i]) / 4.0) / m_SquareGradientDelta;
        }
    }

    template <typename T1, typename T2>
    void
    InitializeRightSideCoefs(T1 &outL, T1 &outM, T1 &outR,
                             const T2 &delta, const T1 &inPrevM, const T1 &inPrevR,
                             const T1 &diffsPrev, const T1 &diffs)
    {
        outL = 1.0 / inPrevM;
        outL *= - delta * (5 * diffs - diffsPrev) / (4.0 * m_SquareGradientDelta);

        outM = 1 + 2 * delta * diffs / m_SquareGradientDelta;
        outM -= outL * inPrevR;

        outR = 0;
    }

    template <typename T1, typename T2, unsigned int TDimension>
    void
    InitializeRightSideCoefs(itk::Vector <T1,TDimension> &outL,
                             itk::Vector <T1,TDimension> &outM,
                             itk::Vector <T1,TDimension> &outR,
                             const T2 &delta, const itk::Vector <T1,TDimension> &inPrevM,
                             const itk::Vector <T1,TDimension> &inPrevR,
                             const itk::Vector <T1,TDimension> &diffsPrev,
                             const itk::Vector <T1,TDimension> &diffs)
    {
        for (unsigned int i = 0;i < TDimension;++i)
        {
            outL[i] = 1.0 / inPrevM[i];
            outL[i] *= - delta * (5 * diffs[i] - diffsPrev[i]) / (4.0 * m_SquareGradientDelta);

            outM[i] = 1 + 2 * delta * diffs[i] / m_SquareGradientDelta;
            outM[i] -= outL[i] * inPrevR[i];

            outR[i] = 0;
        }
    }

    template <typename T1>
    void
    ComputeForwardBacwardSubstitution(std::vector <T1> &outs,
                                      std::vector <T1> &scratch, const std::vector <T1> &data,
                                      const std::vector <T1> &l_coefs, const std::vector <T1> &m_coefs,
                                      const std::vector <T1> &r_coefs)
    {
        unsigned int ln = data.size();

        // compute forward substitution inside scratch
        scratch[0] = data[0];
        for (unsigned int i = 1;i < ln;++i)
            scratch[i] = data[i] - l_coefs[i] * scratch[i-1];

        // And compute backward substitution inside outs
        outs[ln-1] = scratch[ln-1] / m_coefs[ln-1];

        for (int i = ln-2;i >= 0;i--)
            outs[i] = (scratch[i] - r_coefs[i] * outs[i+1]) / m_coefs[i];
    }


    template <typename T1, unsigned int TDimension>
    void
    ComputeForwardBacwardSubstitution(std::vector <itk::Vector <T1,TDimension> > &outs,
                                      std::vector <itk::Vector <T1,TDimension> > &scratch,
                                      const std::vector <itk::Vector <T1,TDimension> > &data,
                                      const std::vector <itk::Vector <T1,TDimension> > &l_coefs,
                                      const std::vector <itk::Vector <T1,TDimension> > &m_coefs,
                                      const std::vector <itk::Vector <T1,TDimension> > &r_coefs)
    {
        unsigned int ln = data.size();

        // compute forward substitution inside scratch
        for (unsigned int j = 0;j < TDimension;++j)
            scratch[0][j] = data[0][j];

        for (unsigned int i = 1;i < ln;++i)
        {
            for (unsigned int j = 0;j < TDimension;++j)
                scratch[i][j] = data[i][j] - l_coefs[i][j] * scratch[i-1][j];
        }

        // And compute backward substitution inside outs
        for (unsigned int j = 0;j < TDimension;++j)
            outs[ln-1][j] = scratch[ln-1][j] / m_coefs[ln-1][j];

        for (int i = ln-2;i >= 0;i--)
        {
            for (unsigned int j = 0;j < TDimension;++j)
                outs[i][j] = (scratch[i][j] - r_coefs[i][j] * outs[i+1][j]) / m_coefs[i][j];
        }
    }

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(InhomogeneousAOSLineDiffusionImageFilter);

    /** Direction in which the filter is to be applied
         * this should be in the range [0,ImageDimension-1]. */
    unsigned int m_Direction;

    DiffusionScalarsImagePointer m_DiffusionScalarsImage;
    double m_TimeStep;
    double m_SquareGradientDelta;
};

} // end namespace anima

#include "animaInhomogeneousAOSLineDiffusionImageFilter.hxx"
