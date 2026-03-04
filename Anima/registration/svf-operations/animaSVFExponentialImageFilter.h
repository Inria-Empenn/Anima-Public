#pragma once

#include <itkImageToImageFilter.h>

namespace anima {

/**
 * @brief Computes the exponentiation of a stationary velocity field using
 * sclaing and squaring and approximated exponential integrators
 *
 * Depending on the order set (0 or 1), the approximated exponentiation for a
 * small enough field is performed using Ferraris et al. ss_aei or Arsigny et
 * al. original approach
 *
 * S. Ferraris et al. Accurate small deformation exponential approximant to
 * integrate large velocity fields: Application to image registration. WBIR 2016
 * V. Arsigny et al. A Log-Euclidean Framework for Statistics on
 * Diffeomorphisms. MICCAI 2006.
 */
template <typename TPixelType, unsigned int Dimension>
class SVFExponentialImageFilter
    : public itk::ImageToImageFilter<
          itk::Image<itk::Vector<TPixelType, Dimension>, Dimension>,
          itk::Image<itk::Vector<TPixelType, Dimension>, Dimension>> {
public:
  using Self = SVFExponentialImageFilter;
  using InputImageType =
      itk::Image<itk::Vector<TPixelType, Dimension>, Dimension>;
  using JacobianImageType =
      itk::Image<itk::Vector<TPixelType, Dimension * Dimension>, Dimension>;
  using OutputImageType =
      itk::Image<itk::Vector<TPixelType, Dimension>, Dimension>;
  using Superclass = itk::ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  itkTypeMacro(SVFExponentialImageFilter, itk::ImageToImageFilter);

  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using JacobianImagePointer = typename JacobianImageType::Pointer;
  using JacobianPixelType = typename JacobianImageType::PixelType;

  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(ExponentiationOrder, unsigned int);
  itkSetMacro(MaximalDisplacementAmplitude, double);

protected:
  SVFExponentialImageFilter() {
    m_ExponentiationOrder = 0;
    m_MaximalDisplacementAmplitude = 0.25;
    m_FieldJacobian = 0;
  }

  virtual ~SVFExponentialImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
  void AfterThreadedGenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(SVFExponentialImageFilter);

  //! Maximal voxel displacement amplitude to be achieved after scaling
  double m_MaximalDisplacementAmplitude;

  //! Exponentiation order (0: Arsigny et al., 1: ss_aei from Ferraris et al.)
  double m_ExponentiationOrder;

  //! Jacobian field (computed only if order 1)
  JacobianImagePointer m_FieldJacobian;

  //! Internal variable that holds the automatically computed number of
  //! recursive squarings
  unsigned int m_NumberOfSquarings;
};

} // end namespace anima

#include "animaSVFExponentialImageFilter.hxx"
