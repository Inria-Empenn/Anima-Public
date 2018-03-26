namespace anima
{

template <class ImageType, class DataImageType>
NonLocalT2DistributionPatchSearcher <ImageType, DataImageType>
::NonLocalT2DistributionPatchSearcher()
{
    m_BetaParameter = 1.0;

    m_MeanMinThreshold = 0.95;
    m_VarMinThreshold = 0.5;
}

template <class ImageType, class DataImageType>
bool
NonLocalT2DistributionPatchSearcher <ImageType, DataImageType>
::TestPatchConformity(unsigned int index, const IndexType &refIndex, const IndexType &movingIndex)
{
    unsigned int numTestImages = m_MeanImages.size();

    for (unsigned int i = 0;i < numTestImages;++i)
    {
        double refMeanValue = m_MeanImages[i]->GetPixel(refIndex);
        double floMeanValue = m_MeanImages[i]->GetPixel(movingIndex);

        double refVarValue = m_VarImages[i]->GetPixel(refIndex);
        double floVarValue = m_VarImages[i]->GetPixel(movingIndex);

        double meanRate = refMeanValue / floMeanValue;
        double varianceRate = refVarValue / floVarValue;

        bool testRate = ( meanRate > m_MeanMinThreshold ) && ( meanRate < ( 1.0 / m_MeanMinThreshold ) ) &&
                ( varianceRate > m_VarMinThreshold ) && ( varianceRate < ( 1.0 / m_VarMinThreshold ) );

        if (!testRate)
            return false;
    }

    return true;
}

template <class ImageType, class DataImageType>
double
NonLocalT2DistributionPatchSearcher <ImageType, DataImageType>
::ComputeWeightValue(unsigned int index, ImageRegionType &refPatch, ImageRegionType &movingPatch)
{
    typedef itk::ImageRegionConstIteratorWithIndex< DataImageType > InIteratorType;

    unsigned int numTestImages = m_PatchTestImages.size();
    double globalWeightValue = 0.0;

    for (unsigned int i = 0;i < numTestImages;++i)
    {
        InIteratorType tmpIt (m_PatchTestImages[i], refPatch);
        InIteratorType tmpMovingIt (m_PatchTestImages[i], movingPatch);

        double tmpDiffValue;

        double weightValue = 0.0;
        unsigned int numVoxels = 0;

        while (!tmpIt.IsAtEnd())
        {
            tmpDiffValue = (double)tmpIt.Get() - (double)tmpMovingIt.Get();
            weightValue += tmpDiffValue * tmpDiffValue;

            ++numVoxels;
            ++tmpIt;
            ++tmpMovingIt;
        }

        globalWeightValue -= weightValue / (2.0 * m_BetaParameter * m_NoiseCovariances[i] * numVoxels);
    }

    return std::exp(globalWeightValue / numTestImages);
}

} // end namespace anima
