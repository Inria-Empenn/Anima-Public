#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkImage.h>

namespace anima
{

/**
 * @brief Compute the Jacobian matrix in real coordinates of a displacement field
 *
 * The Jacobian matrix is computed as a linear least squares problem based on multiple directional derivatives around each pixel
 * The output vector in each voxel is of size Dimension^2 and is stored on rows first, i.e.
 * J[0] = dT_0/dx_0, J[1] = dT_0/dx_1, J[2] = dT_0/dx_2, J[3] = dT_1/dx_0, ...
 */
template <typename TPixelType, typename TOutputPixelType, unsigned int Dimension>
class JacobianMatrixImageFilter :
public itk::ImageToImageFilter< itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> ,
        itk::Image <itk::Vector <TOutputPixelType, Dimension * Dimension>, Dimension> >
{
public:
    typedef JacobianMatrixImageFilter Self;
    typedef typename itk::Image <itk::Vector <TPixelType, Dimension>, Dimension> InputImageType;
    typedef typename itk::Image <itk::Vector <TOutputPixelType, Dimension * Dimension>, Dimension> OutputImageType;
    typedef typename itk::Image <TOutputPixelType, Dimension> DeterminantImageType;
    typedef itk::ImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(JacobianMatrixImageFilter, itk::ImageToImageFilter)

    typedef typename InputImageType::IndexType IndexType;
    typedef typename InputImageType::RegionType RegionType;
    typedef typename InputImageType::PointType PointType;
    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction <InputImageType> InterpolatorType;
    typedef typename InterpolatorType::Pointer InterpolatorPointer;

    itkSetMacro(NoIdentity, bool)
    itkSetMacro(ComputeDeterminant, bool)
    itkGetObjectMacro(DeterminantImage, DeterminantImageType)

    void SetNeighborhood(unsigned int val);

protected:
    JacobianMatrixImageFilter()
    {
        m_NoIdentity = false;
        m_Neighborhood = 1;
        m_OnlySixConnectivity = true;

        m_ComputeDeterminant = false;
        m_DeterminantImage = 0;
    }

    virtual ~JacobianMatrixImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    bool CheckFaceConnectivity(const IndexType &internalIndex, const IndexType &currentIndex);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(JacobianMatrixImageFilter);

    //! Add identity to the obtained Jacobian matrix
    bool m_NoIdentity;

    //! Neighborhood from which to compute Jacobian matrix, 0 -> do only 6 connectivity approximation
    unsigned int m_Neighborhood;
    bool m_OnlySixConnectivity;

    typename DeterminantImageType::Pointer m_DeterminantImage;
    bool m_ComputeDeterminant;

    InterpolatorPointer m_FieldInterpolator;
};

} // end namespace anima

#include "animaJacobianMatrixImageFilter.hxx"
