#pragma once

#include <itkInterpolateImageFunction.h>
#include <animaMultiCompartmentModel.h>
#include <mutex>
#include <animaMCMWeightedAverager.h>

namespace anima
{

template <class TInputImage, class TCoordRep = double>
class MCMLinearInterpolateImageFunction :
public itk::InterpolateImageFunction<TInputImage,TCoordRep>
{
public:
    /** Standard class typedefs. */
    typedef MCMLinearInterpolateImageFunction Self;
    typedef itk::InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MCMLinearInterpolateImageFunction,
                 InterpolateImageFunction)

    /** InputImageType typedef support. */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename TInputImage::PixelType PixelType;
    typedef typename Superclass::RealType RealType;
    typedef typename Superclass::SizeType SizeType;

    /** Dimension underlying input image. */
    itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

    /** Index typedef support. */
    typedef typename Superclass::IndexType       IndexType;
    typedef typename Superclass::IndexValueType  IndexValueType;

    /** ContinuousIndex typedef support. */
    typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

    typedef typename Superclass::OutputType OutputType;

    // Multi-compartment models typedefs
    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    typedef anima::MCMWeightedAverager AveragerType;
    typedef AveragerType::Pointer AveragerPointer;

    /** Evaluate the function at a ContinuousIndex position
     *
     * Returns the linearly interpolated image intensity at a
     * specified point position. No bounds checking is done.
     * The point is assumed to lie within the image buffer.
     *
     * ImageFunction::IsInsideBuffer() can be used to check bounds before
     * calling the method. */
    virtual OutputType EvaluateAtContinuousIndex(const ContinuousIndexType & index ) const ITK_OVERRIDE;

    void SetReferenceOutputModel(MCModelPointer &model);

    virtual void SetInputImage(const InputImageType *ptr) ITK_OVERRIDE;

    std::vector <AveragerPointer> &GetAveragers() const {return m_MCMAveragers;}

    SizeType GetRadius() const override
    {
        return SizeType::Filled(1);
    }

protected:
    MCMLinearInterpolateImageFunction();
    virtual ~MCMLinearInterpolateImageFunction() {}

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

    //! Tests if input and output models of the interpolator are compatible
    void TestModelsAdequation(MCModelPointer &inputModel, MCModelPointer &outputModel);

    //! Resets averager pointers, can be overloaded to handle more models
    virtual void ResetAveragePointers(MCModelPointer &model);

    //! Check if model can actually be interpolated
    virtual bool CheckModelCompatibility(MCModelPointer &model);

    //! Sets averager specific parameters if sub-classes are derived
    virtual void SetSpecificAveragerParameters(unsigned int threadIndex) const {}

    unsigned int GetFreeWorkIndex() const;
    void UnlockWorkIndex(unsigned int index) const;

private:
    MCMLinearInterpolateImageFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /** Number of neighbors used in the interpolation */
    static const unsigned long m_Neighbors;
    static const unsigned int m_SphereDimension = 3;

    mutable std::mutex m_LockUsedModels;
    mutable std::vector <int> m_UsedModels;

    mutable std::vector < std::vector <MCModelPointer> > m_ReferenceInputModels;
    mutable std::vector < std::vector <double> > m_ReferenceInputWeights;
    mutable std::vector <AveragerPointer> m_MCMAveragers;
};

} // end namespace anima

#include "animaMCMLinearInterpolateImageFunction.hxx"
