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
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef itk::Image <PixelType, NDimensions> ScalarImageType;
    typedef itk::Image <double, NDimensions> WeightImageType;
    typedef typename WeightImageType::Pointer WeightImagePointer;
    typedef typename ScalarImageType::IndexType IndexType;
    typedef typename ScalarImageType::PointType PointType;
    typedef itk::VectorImage <PixelType, NDimensions> VectorImageType;
    typedef itk::Image <unsigned char, NDimensions> MaskImageType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(BlockMatchingInitializer, itk::Object)

    typedef typename ScalarImageType::RegionType ImageRegionType;

    typedef typename ScalarImageType::Pointer ScalarImagePointer;
    typedef typename VectorImageType::Pointer VectorImagePointer;
    typedef typename MaskImageType::Pointer MaskImagePointer;

    itkSetMacro(NumberOfThreads, unsigned int)
    itkGetMacro(NumberOfThreads, unsigned int)

    void SetPercentageKept(double val);
    itkGetMacro(PercentageKept, double)

    void SetScalarVarianceThreshold(double val);
    itkGetMacro(ScalarVarianceThreshold, double)

    void SetOrientedModelVarianceThreshold(double val);
    itkGetMacro(OrientedModelVarianceThreshold, double)

    void SetRequestedRegion(const ImageRegionType &val);
    itkGetMacro(RequestedRegion, ImageRegionType)

    void SetBlockSpacing(unsigned int val);
    void SetBlockSize(unsigned int val);

    itkGetMacro(BlockSpacing, unsigned int)
    itkGetMacro(BlockSize, unsigned int)

    void clearGenerationMasks() {m_GenerationMasks.clear();}
    void AddGenerationMask(MaskImageType *mask);

    virtual void AddReferenceImage(itk::ImageBase <NDimensions> *refImage);
    itk::ImageBase <NDimensions> *GetFirstReferenceImage();

    ScalarImageType *GetReferenceScalarImage(unsigned int i) {return m_ReferenceScalarImages[i];}
    VectorImageType *GetReferenceVectorImage(unsigned int i) {return m_ReferenceVectorImages[i];}

    void Update();

    // Does the splitting and calls RegionBlockGenerator on a sub region
    static ITK_THREAD_RETURN_TYPE ThreadBlockGenerator(void *arg);

    std::vector <ImageRegionType> &GetOutput();
    std::vector <unsigned int> &GetMaskStartingIndexes();

    std::vector <PointType> &GetOutputPositions();

    struct BlockGeneratorThreadStruct
    {
        Pointer Filter;
        std::vector < std::vector <ImageRegionType> > tmpOutput;
        std::vector <unsigned int> blockStartOffsets;
        std::vector < std::vector <unsigned int> > startBlocks, nb_blocks;
        std::vector <unsigned int> totalNumberOfBlocks;
        std::vector < std::vector <double> > blocks_variances;
        unsigned int maskIndex;
        std::vector < std::vector <PointType> > blocks_positions;

        //! To be able to inherit from it
        virtual ~BlockGeneratorThreadStruct() {}
    };

    void RegionBlockGenerator(BlockGeneratorThreadStruct *workStr, unsigned int threadId);

protected:
    BlockMatchingInitializer() : Superclass()
    {
        m_BlockSize = 5;
        m_BlockSpacing = 3;
        m_NumberOfThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();

        m_ScalarVarianceThreshold = 5.0;
        m_OrientedModelVarianceThreshold = 0.0;

        m_PercentageKept = 0.8;
        m_GenerationMasks.clear();

        m_UpToDate = false;
    }

    virtual ~BlockMatchingInitializer() {}

    void ComputeBlocksOnGenerationMask(unsigned int maskIndex);

    virtual void InitializeThreading(unsigned int maskIndex, BlockGeneratorThreadStruct *&workStr);

    bool CheckBlockConditions(ImageRegionType &region, double &blockVariance, BlockGeneratorThreadStruct *workStr,
                              unsigned int threadId);

    bool CheckScalarVariance(ScalarImageType *refImage, ImageRegionType &region, double &blockVariance);
    virtual bool CheckOrientedModelVariance(unsigned int imageIndex, ImageRegionType &region, double &blockVariance,
                                            BlockGeneratorThreadStruct *workStr, unsigned int threadId);

    bool ProgressCounter(std::vector <unsigned int> &counter, std::vector <unsigned int> &bounds);

    struct pair_comparator
    {
        bool operator() (const std::pair<double, std::pair <PointType, ImageRegionType> > & f,
                         const std::pair<double, std::pair <PointType, ImageRegionType> > & s)
        { return (f.first < s.first); }
    };

private:
    BlockMatchingInitializer(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    unsigned int m_NumberOfThreads;

    double m_ScalarVarianceThreshold;
    double m_OrientedModelVarianceThreshold;
    double m_PercentageKept;

    ImageRegionType m_RequestedRegion;

    std::vector <ScalarImagePointer> m_ReferenceScalarImages;
    std::vector <VectorImagePointer> m_ReferenceVectorImages;

    std::vector <MaskImagePointer> m_GenerationMasks;

    std::vector <ImageRegionType> m_Output;
    std::vector <unsigned int> m_MaskStartingIndexes;
    std::vector <PointType> m_OutputPositions;

    bool m_UpToDate;
};

} // end of namespace anima

#include "animaBlockMatchInitializer.hxx"
