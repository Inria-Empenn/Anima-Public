#pragma once

#include <itkImageRegion.h>
#include <itkObject.h>
#include <itkMultiThreader.h>

#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{

template <class PixelType, unsigned int NDimensions=3>
class BlockMatchingInitializer
        : public itk::Object
{
public:
    typedef itk::Object Superclass;
    typedef BlockMatchingInitializer <PixelType,NDimensions> Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef itk::Image <PixelType, NDimensions> ScalarImageType;
    typedef typename ScalarImageType::IndexType IndexType;
    typedef itk::VectorImage <PixelType, NDimensions> VectorImageType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(BlockMatchingInitializer, Object);

    typedef typename ScalarImageType::RegionType ImageRegionType;

    typedef typename ScalarImageType::Pointer ScalarImagePointer;
    typedef typename VectorImageType::Pointer VectorImagePointer;

    itkSetMacro(NumberOfThreads, unsigned int);
    itkGetMacro(NumberOfThreads, unsigned int);

    itkSetMacro(PercentageKept, double);
    itkGetMacro(PercentageKept, double);

    itkSetMacro(ScalarVarianceThreshold, double);
    itkGetMacro(ScalarVarianceThreshold, double);

    itkSetMacro(TensorVarianceThreshold, double);
    itkGetMacro(TensorVarianceThreshold, double);

    itkSetMacro(RequestedRegion, ImageRegionType);
    itkGetMacro(RequestedRegion, ImageRegionType);

    itkSetMacro(BlockSpacing, unsigned int);
    itkSetMacro(BlockSize, unsigned int);

    itkGetMacro(BlockSpacing, unsigned int);
    itkGetMacro(BlockSize, unsigned int);

    itkSetMacro(ComputeOuterDam, bool);
    itkGetMacro(ComputeOuterDam, bool);

    itkSetMacro(DamDistance, double);
    itkGetMacro(DamDistance, double);

    std::vector <IndexType> &GetDamIndexes() {return m_DamIndexes;}

    void AddReferenceImage(itk::ImageBase <NDimensions> *refImage);
    itk::ImageBase <NDimensions> *GetFirstReferenceImage();

    void Update();

    // Does the splitting and calls RegionBlockGenerator on a sub region
    static ITK_THREAD_RETURN_TYPE ThreadBlockGenerator(void *arg);

    void RegionBlockGenerator(std::vector <unsigned int> &num_start, std::vector <unsigned int> &num_blocks,
                              std::vector <ImageRegionType> &tmpOutput, std::vector <double> &block_variances,
                              unsigned int &total_num_blocks);

    std::vector <ImageRegionType> &GetOutput();

protected:
    BlockMatchingInitializer() : Superclass()
    {
        m_BlockSize = 5;
        m_BlockSpacing = 3;
        m_NumberOfThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();

        m_ScalarVarianceThreshold = 5.0;
        m_TensorVarianceThreshold = 0.0;

        m_PercentageKept = 0.8;

        m_ComputeOuterDam = false;
        m_DamDistance = 4;
        m_DamIndexes.clear();
    }

    virtual ~BlockMatchingInitializer() {}

    void ComputeOuterDamFromBlocks();

    bool CheckBlockConditions(ImageRegionType &region, double &blockVariance);
    bool CheckScalarVariance(ScalarImageType *refImage, ImageRegionType &region, double &blockVariance);
    bool CheckTensorVariance(VectorImageType *refImage, ImageRegionType &region, double &blockVariance);

    struct BlockGeneratorThreadStruct
    {
        Pointer Filter;
        std::vector < std::vector <ImageRegionType> > tmpOutput;
        std::vector < std::vector <unsigned int> > startBlocks, nb_blocks;
        std::vector <unsigned int> totalNumberOfBlocks;
        std::vector < std::vector <double> > blocks_variances;
    };

    struct pair_comparator
    {
        bool operator() (const std::pair<double, ImageRegionType> & f, const std::pair<double, ImageRegionType> & s)
        { return (f.first < s.first); }
    };

private:
    BlockMatchingInitializer(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    unsigned int m_NumberOfThreads;

    double m_ScalarVarianceThreshold;
    double m_TensorVarianceThreshold;
    double m_PercentageKept;

    bool m_ComputeOuterDam;
    //! Distance from blocks to the dam, is in voxel coordinates !
    double m_DamDistance;
    std::vector <IndexType> m_DamIndexes;

    ImageRegionType m_RequestedRegion;

    std::vector <ScalarImagePointer> m_ReferenceScalarImages;
    std::vector <VectorImagePointer> m_ReferenceVectorImages;

    std::vector <ImageRegionType> m_Output;
};

} // end of namespace anima

#include "animaBlockMatchInitializer.hxx"
