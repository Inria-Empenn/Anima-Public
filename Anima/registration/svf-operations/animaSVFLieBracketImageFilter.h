#pragma once

#include <itkImageToImageFilter.h>

namespace anima {

/**
 * @brief Computes the Lie bracket between two fields u and v as expressed by
 * Bossa et al.
 *
 * There is a discrepancy in between Vercauteren et al. and Bossa et al.
 * formulation. This Lie bracket implements the Vercauteren et al. formulation
 * [u,v](x) = Jac(u)(x).v(x) - Jac(v)(x).u(x)
 * M. Bossa et al. "Contributions to 3D diffeomorphic atlas estimation :
 * application to brain images.", MICCAI 2007, p. 667–674. T. Vercauteren et al.
 * "Symmetric Log-Domain Diffeomorphic Registration: A Demons-based Approach.",
 * MICCAI 2008, p. 754-761.
 */
template <typename TPixelType, unsigned int Dimension>
class SVFLieBracketImageFilter
    : public itk::ImageToImageFilter<
          itk::Image<itk::Vector<TPixelType, Dimension>, Dimension>,
          itk::Image<itk::Vector<TPixelType, Dimension>, Dimension>> {
public:
  using Self = SVFLieBracketImageFilter;
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

  itkTypeMacro(SVFLieBracketImageFilter, itk::ImageToImageFilter);

  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using JacobianImagePointer = typename JacobianImageType::Pointer;
  using JacobianPixelType = typename JacobianImageType::PixelType;

  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetObjectMacro(FirstFieldJacobian, JacobianImageType);
  itkSetObjectMacro(SecondFieldJacobian, JacobianImageType);
  itkGetObjectMacro(FirstFieldJacobian, JacobianImageType);
  itkGetObjectMacro(SecondFieldJacobian, JacobianImageType);

protected:
  SVFLieBracketImageFilter() {
    this->SetNumberOfRequiredInputs(2);
    m_FirstFieldJacobian = 0;
    m_SecondFieldJacobian = 0;
  }

  virtual ~SVFLieBracketImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(SVFLieBracketImageFilter);

  JacobianImagePointer m_FirstFieldJacobian, m_SecondFieldJacobian;
};

} // end namespace anima

#include "animaSVFLieBracketImageFilter.hxx"
