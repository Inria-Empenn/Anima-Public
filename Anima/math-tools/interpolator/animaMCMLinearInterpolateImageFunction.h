#pragma once

#include <animaMCMWeightedAverager.h>
#include <animaMultiCompartmentModel.h>
#include <itkInterpolateImageFunction.h>
#include <mutex>

namespace anima {

template <class TInputImage, class TCoordRep = double>
class MCMLinearInterpolateImageFunction
    : public itk::InterpolateImageFunction<TInputImage, TCoordRep> {
public:
  /** Standard class typedefs. */
  using Self = MCMLinearInterpolateImageFunction;
  using Superclass = itk::InterpolateImageFunction<TInputImage, TCoordRep>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MCMLinearInterpolateImageFunction, InterpolateImageFunction);

  /** InputImageType typedef support. */
  using InputImageType = typename Superclass::InputImageType;
  using PixelType = typename TInputImage::PixelType;
  using RealType = typename Superclass::RealType;
  using SizeType = typename Superclass::SizeType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Index typedef support. */
  using IndexType = typename Superclass::IndexType;
  using IndexValueType = typename Superclass::IndexValueType;

  /** ContinuousIndex typedef support. */
  using ContinuousIndexType = typename Superclass::ContinuousIndexType;

  using OutputType = typename Superclass::OutputType;

  // Multi-compartment models typedefs
  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = MCModelType::Pointer;

  using AveragerType = anima::MCMWeightedAverager;
  using AveragerPointer = AveragerType::Pointer;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assumed to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(
      const ContinuousIndexType &index) const ITK_OVERRIDE;

  void SetReferenceOutputModel(MCModelPointer &model);

  virtual void SetInputImage(const InputImageType *ptr) ITK_OVERRIDE;

  std::vector<AveragerPointer> &GetAveragers() const { return m_MCMAveragers; }

  SizeType GetRadius() const override { return SizeType::Filled(1); }

  itkSetMacro(DDIInterpolationMethod, unsigned int);

protected:
  MCMLinearInterpolateImageFunction();
  virtual ~MCMLinearInterpolateImageFunction() {}

  template <class T>
  bool isZero(const itk::VariableLengthVector<T> &value) const {
    for (unsigned int i = 0; i < value.GetNumberOfElements(); ++i) {
      if (value[i] != 0)
        return false;
    }

    return true;
  }

  // Fake method for compilation purposes, should never go in there
  template <class T> bool isZero(T &data) const {
    itkExceptionMacro("Access to unauthorized method");
    return true;
  }

  //! Tests if input and output models of the interpolator are compatible
  void TestModelsAdequation(MCModelPointer &inputModel,
                            MCModelPointer &outputModel);

  //! Resets averager pointers, can be overloaded to handle more models
  void ResetAveragePointers(MCModelPointer &model);

  //! Check if model can actually be interpolated
  bool CheckModelCompatibility(MCModelPointer &model);

  //! Sets averager specific parameters
  void SetSpecificAveragerParameters(unsigned int threadIndex) const;

  unsigned int GetFreeWorkIndex() const;
  void UnlockWorkIndex(unsigned int index) const;

private:
  MCMLinearInterpolateImageFunction(const Self &); // purposely not implemented
  void operator=(const Self &);                    // purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long m_Neighbors;
  static const unsigned int m_SphereDimension = 3;

  unsigned int m_DDIInterpolationMethod;

  mutable std::mutex m_LockUsedModels;
  mutable std::vector<int> m_UsedModels;

  mutable std::vector<std::vector<MCModelPointer>> m_ReferenceInputModels;
  mutable std::vector<std::vector<double>> m_ReferenceInputWeights;
  mutable std::vector<AveragerPointer> m_MCMAveragers;
};

} // end namespace anima

#include "animaMCMLinearInterpolateImageFunction.hxx"
