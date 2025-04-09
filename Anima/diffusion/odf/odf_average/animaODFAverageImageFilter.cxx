#include "animaODFAverageImageFilter.h"
#include <animaLogExpMapsUnitSphere.h>
#include <animaVectorOperations.h>

#include <itkImageRegionIterator.h>

#define NB_SAMPLES_THETA 10
#define NB_SAMPLES_PHI   (2 * NB_SAMPLES_THETA)
#define DELTA_PHI        (2.0 * M_PI / static_cast<double>(NB_SAMPLES_PHI))
#define DELTA_THETA      (M_PI / (static_cast<double>(NB_SAMPLES_THETA) - 1.0))


namespace anima
{
    ODFAverageImageFilter::ODFAverageImageFilter()
    {
        m_EpsValueTestNormalized = std::pow(10, -6);
        m_EpsValueTestNull = std::pow(10, -12);
        m_VectorLength = 0;
        m_ODFSHBasis = nullptr;
        //m_SamplePoints = this->GetSamplePoints();
    }

    ODFAverageImageFilter::~ODFAverageImageFilter() = default;

    void ODFAverageImageFilter::AddWeightImage(const unsigned int i, const WeightImagePointer &weightImage)
    {
        /**
         * Add a newly read weight image to the vector m_WeightImages 
         */
        if (i == m_WeightImages.size())
        {
            m_WeightImages.push_back(weightImage);
            return;
        }

        if (i > m_WeightImages.size())
            itkExceptionMacro("Weight images must be added contiguously.");

        m_WeightImages[i] = weightImage;
    }


    void ODFAverageImageFilter::BeforeThreadedGenerateData()
    {
        /**
         * Preliminary computations
         */
        Superclass::BeforeThreadedGenerateData();

        m_VectorLength = this->GetInput(0)->GetVectorLength();
        //Get the SH order from the vector length
        unsigned int odfSHOrder = std::round(-1.5 + 0.5 * std::sqrt(8.0 * static_cast<double>(m_VectorLength) + 1.0));

        //Discretize SH (m_SpherHarm will contain the values of each spherical harmonic (columns) for each sample point (lines))
        m_SpherHarm.set_size(NB_SAMPLES_PHI * NB_SAMPLES_THETA, m_VectorLength);
        m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(odfSHOrder);

        double sqrt2 = std::sqrt(2.0);
        unsigned int k = 0;
        for (unsigned int i = 0; i < NB_SAMPLES_THETA; ++i)
        {
            double theta = static_cast<double>(i) * DELTA_THETA;

            for (unsigned int j = 0; j < NB_SAMPLES_PHI; ++j)
            {
                double phi = static_cast<double>(j) * DELTA_PHI;
                unsigned int c = 0;
                for (double l = 0; l <= odfSHOrder; l += 2)
                    for (double m = -l; m <= l; m++)
                        m_SpherHarm.put(k, c++, m_ODFSHBasis->getNthSHValueAtPosition(l, m, theta, phi));
                k++;
            }
        }

        m_SolveSHMatrix = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose();
    }


    void ODFAverageImageFilter::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
    {
        /**
         *Main function. For each voxel of output image, compute the average of the corresponding input voxels in input images. 
         */

        //Initializations
        using InputImageIteratorType = itk::ImageRegionConstIterator<InputImageType>;
        using InputWeightIteratorType = itk::ImageRegionConstIterator<WeightImageType>;
        using OutputImageIteratorType = itk::ImageRegionIterator<OutputImageType>;

        unsigned int numInputs = this->GetNumberOfIndexedInputs();

        std::vector<InputImageIteratorType> inItrs(numInputs); //vector of input images (which are here iterators on the input images)
        std::vector<InputWeightIteratorType> weightItrs(numInputs); //vector of iterators on weight images
        for (unsigned int i = 0; i < numInputs; ++i)
        {
            inItrs[i] = InputImageIteratorType(this->GetInput(i), outputRegionForThread);
            weightItrs[i] = InputWeightIteratorType(m_WeightImages[i], outputRegionForThread);
        }

        OutputImageIteratorType outItr(this->GetOutput(), outputRegionForThread); //iterator on output image

        VectorType coefODFValue(m_VectorLength);
        VectorType sampledODFValue(NB_SAMPLES_PHI * NB_SAMPLES_THETA);
        VectorArrayType currentODFs(numInputs);
        VectorArrayType notNullCurrentODFs;
        VectorType currentWeights(numInputs);
        VectorType currentWeightsNotNullODFs;
        VectorType inputValue(m_VectorLength);
        OutputPixelType outputValue(m_VectorLength);

        //Average computation
        while (!outItr.IsAtEnd())
        {
            double weightSum = 0.0;
            double normValue = 1.0;
            int numNotNullODFs;
            int indexODFtoNormalize;
            for (unsigned int i = 0; i < numInputs; ++i)
            {
                double weightValue = weightItrs[i].Get();
                weightSum += weightValue;
                currentWeights[i] = weightValue;
            }

            if (weightSum < m_EpsValueTestNull)
            //If in a voxel, the sum of weights is equal (or very close) to 0, then we output a null distribution
            {
                outputValue.Fill(0.0);
                outItr.Set(outputValue);

                //Move to next voxel
                for (unsigned int i = 0; i < numInputs; ++i)
                {
                    ++inItrs[i];
                    ++weightItrs[i];
                }

                ++outItr;
                this->IncrementNumberOfProcessedPoints();
            }

            else
            {
                //std::cout<<"1"<<std::endl;
                for (unsigned int i = 0; i < numInputs; ++i)
                {
                    this->ConvertODF(inItrs[i].Get(), inputValue);
                    currentODFs[i] = inputValue;
                    //std::cout<<"2"<<std::endl;

                }
                numNotNullODFs = this->RemoveNullODFs(currentODFs, currentWeights, numInputs, notNullCurrentODFs, currentWeightsNotNullODFs);
                //std::cout<<"3"<<std::endl;
                switch (numNotNullODFs)
                {
                    case 0:
                        //std::cout<<"case 0"<<std::endl;
                        outputValue.Fill(0.0);
                        break;
                    
                    case 1:
                        //std::cout<<"case 1"<<std::endl;
                        for (unsigned int i = 0; i < m_VectorLength; i++){
                            outputValue[i] = notNullCurrentODFs[0][i];
                        }
                        break;
                    
                    default:
                        //std::cout<<"default case"<<std::endl;
                        indexODFtoNormalize = this->GetIndexODFtoNormalize(notNullCurrentODFs, numNotNullODFs);
                        //std::cout<<indexODFtoNormalize<<std::endl;
                        for (unsigned int i = 0; i < numNotNullODFs; i++)
                        {
                            /*
                            std::cout<<"Input ODF "<<i<<std::endl;
                            std::cout<<"[";
                            for (double v : notNullCurrentODFs[i])
                            {
                                std::cout<<v<<", ";
                            }
                            std::cout<<"]"<<std::endl;
                            */
                            this->DiscretizeODF(notNullCurrentODFs[i], sampledODFValue);
                            coefODFValue = this->GetSquareRootODFCoef(sampledODFValue);
                            //std::cout<<"5"<<std::endl;
                            if (i == indexODFtoNormalize)
                            {
                                /*
                                std::cout<<"Before normalization "<<i<<std::endl;
                                std::cout<<"[";
                                for (double v : coefODFValue)
                                {
                                    std::cout<<v<<", ";
                                }
                                std::cout<<"]"<<std::endl;
                                */
                                normValue = this->NormalizeODF(coefODFValue, notNullCurrentODFs[i]);
                                /*
                                std::cout<<"Normalized ODF "<<i<<std::endl;
                                std::cout<<"[";
                                for (double v : notNullCurrentODFs[i])
                                {
                                    std::cout<<v<<", ";
                                }
                                std::cout<<"]"<<std::endl;
                                std::cout<<"index "<<i<<std::endl;
                                std::cout<<normValue<<std::endl;
                                */
                            }
                            else
                            {
                                notNullCurrentODFs[i] = coefODFValue;
                                //std::cout<<"not index "<<i<<std::endl;
                            }
                            //std::cout<<std::endl;
                            //std::cout<<std::endl;
                        }
                        coefODFValue = this->GetAverageHisto(notNullCurrentODFs, currentWeightsNotNullODFs, numNotNullODFs);
                        /*
                        std::cout<<"Average ODF"<<std::endl;
                        std::cout<<"[";
                        for (double v : coefODFValue)
                        {
                            std::cout<<v<<", ";
                        }
                        std::cout<<"]"<<std::endl;
                        */
                        //std::cout<<"6"<<std::endl;
                        if (indexODFtoNormalize != -1)
                        {
                            this->InverseNormalization(normValue, coefODFValue);
                            //std::cout<<"7"<<std::endl;
                        }
                        this->DiscretizeODF(coefODFValue, sampledODFValue);
                        outputValue = this->GetSquareODFCoef(sampledODFValue);
                        //std::cout<<"8"<<std::endl;
                }

                notNullCurrentODFs.clear();
                currentWeightsNotNullODFs.clear();

                outItr.Set(outputValue);

                //Move to next voxel
                for (unsigned int i = 0; i < numInputs; ++i)
                {
                    ++inItrs[i];
                    ++weightItrs[i];
                }

                ++outItr;
                this->IncrementNumberOfProcessedPoints();
 
            }
        }
    }
            

    void ODFAverageImageFilter::AfterThreadedGenerateData()
    {
        /**
         * Delete pointers which are not SmartPointers
         */
        delete m_ODFSHBasis;
    }

    void ODFAverageImageFilter::ConvertODF(const InputPixelType &inputODF, VectorType &convertedODF)
    {
        for (unsigned int i = 0; i < m_VectorLength; i++){
            convertedODF[i] = inputODF[i]; 
        }
    }

    int ODFAverageImageFilter::RemoveNullODFs(const VectorArrayType &allODFs, const VectorType &allweights, int nbODFs, VectorArrayType &selectedODFs, VectorType &selectedWeights)
    {
        int numNotNullODFs = 0;
        for (unsigned int i = 0; i < nbODFs; i++){
            VectorType odf = allODFs[i];
            bool currentODFisNull = true;
            int j = 0;
            while (currentODFisNull && j<m_VectorLength)
            {
                if (!std::abs(odf[j])<m_EpsValueTestNull)
                {
                    currentODFisNull = false;
                }
                j++;
            }
            if (!currentODFisNull)
            {
                selectedODFs.push_back(odf);
                selectedWeights.push_back(allweights[i]);
                numNotNullODFs++;
            }
        }
        return numNotNullODFs;
    }

    int ODFAverageImageFilter::GetIndexODFtoNormalize(VectorArrayType &odfs, int &numNotNullODFs)
    {
        int index = -1; //if at the end, index is still -1, it means that not any ODF needs to be normalized
        int i=0;
        while (index == -1 && i < numNotNullODFs)
        {
            if (std::abs(odfs[i][0]-1/std::sqrt(4*M_PI)) >= m_EpsValueTestNormalized)
            {
                index = i;
            }
            i++;
        }
        return index;

    }


    void ODFAverageImageFilter::DiscretizeODF(const VectorType &coefODF, VectorType &sampledODF)
    {
        /**
         * Move an ODF from representation by SH coefficients to a sample form (vector with values at sample points)
         */
        unsigned int k = 0;
        for (unsigned int i = 0; i < NB_SAMPLES_THETA; ++i)
        {
            double theta = static_cast<double>(i) * DELTA_THETA;
            for (unsigned int j = 0; j < NB_SAMPLES_PHI; ++j)
            {
                double phi = static_cast<double>(j) * DELTA_PHI;
                sampledODF[k] = m_ODFSHBasis->getValueAtPosition(coefODF, theta, phi);
                k++;
            }
        }
    }


    double ODFAverageImageFilter::NormalizeODF(const VectorType &inputODF, VectorType &normalizedODF)
    {
        /**
         * Normalizes an ODF so that it integrates to 1 over a sphere, by multiplying each SH component by a constant factor normValue
         */
        double sumSquareCoefs = 0;
        double normValue;
        for (unsigned int i = 0; i < m_VectorLength; i++){
            sumSquareCoefs+=std::pow(inputODF[i],2);
        }
        if (std::abs(sumSquareCoefs)<m_EpsValueTestNull){
            normValue = 0.0;
            for (unsigned int i = 0; i < m_VectorLength; i++){
                normalizedODF[i] = normValue;
            }
        }
        else
        {
            normValue = 1/std::sqrt(sumSquareCoefs);
            for (unsigned int i = 0; i < m_VectorLength; i++){
                normalizedODF[i] = inputODF[i]*normValue;
            }
        }
        return normValue;
    }

    void ODFAverageImageFilter::InverseNormalization(const double &normValue, VectorType &odf)
    {
        if (!std::abs(normValue) < m_EpsValueTestNull)
        {
            for (unsigned int j = 0; j < m_VectorLength; j++){
                odf[j] = odf[j]/normValue;
            }  
        }
    }
    

    ODFAverageImageFilter::VectorType ODFAverageImageFilter::GetSquareRootODFCoef(const VectorType &odf)
    {
        /**
         * Returns the square root of an input ODF
         * Output ODF is represented by coefficients
         * Input ODF is in sample form
         * 
         */

        MatrixType squareRootOdf(NB_SAMPLES_THETA * NB_SAMPLES_PHI, 1);
        MatrixType squareRootCoef(m_VectorLength, 1);

        //Compute square root of each sample point
        for (unsigned int i = 0; i < NB_SAMPLES_THETA * NB_SAMPLES_PHI; ++i)
            squareRootOdf.put(i, 0, std::sqrt(std::max(0.0, odf[i])));

        //Convert from sample form to representation by coefficients
        squareRootCoef = m_SolveSHMatrix * squareRootOdf;

        //Store the result in a vector instead of a matrix
        VectorType squareRootModelValue(m_VectorLength);
        for (unsigned int i = 0; i < m_VectorLength; ++i)
            squareRootModelValue[i] = squareRootCoef.get(i, 0);

        return squareRootModelValue;
    }


    ODFAverageImageFilter::OutputPixelType ODFAverageImageFilter::GetSquareODFCoef(const VectorType &odf)
    {
        /**
         * Returns the square root of an input ODF
         * Output ODF is represented by coefficients
         * Input ODF is in sample form
         * 
         */

        MatrixType squareOdf(NB_SAMPLES_THETA * NB_SAMPLES_PHI, 1);
        MatrixType squareCoef(m_VectorLength, 1);

        //Compute the square of each sample point
        for (unsigned int i = 0; i < NB_SAMPLES_THETA * NB_SAMPLES_PHI; ++i)
            squareOdf.put(i, 0, std::pow(odf[i], 2));

        //Convert from sample form to representation by coefficients
        squareCoef = m_SolveSHMatrix * squareOdf;

        //Store the result in a vector instead of a matrix
        OutputPixelType squareModelValue(m_VectorLength);
        for (unsigned int i = 0; i < m_VectorLength; ++i)
            squareModelValue[i] = squareCoef.get(i, 0);

        return squareModelValue;
    }


    ODFAverageImageFilter::VectorType ODFAverageImageFilter::GetAverageHisto(const VectorArrayType &coefs, const VectorType &weightValues, int nbODFs)
    {
        /**
         * Compute the weighted average of input ODFs using method of Goh et al., 2O11 (Algorithm 1)
         * coefs: Vector of size numInputs, each element being a vector of size m_VectorLength representing an ODF
         * weightValues: Vector of size numInputs with the weights of the input ODFs (each element is a scalar)
         * resCoef: Output ODF represented by coefficients
         */

        VectorType mean, nextMean = coefs[0];
        VectorType tangent(m_VectorLength), workValue(m_VectorLength);

        const unsigned int maxIter = 100;
        const double epsValue = 0.000035;
        unsigned int nIter = 0;
        double normTan = 1.0;

        while (nIter < maxIter && normTan > epsValue)
        {
            mean = nextMean;

            std::fill(tangent.begin(), tangent.end(), 0.0);
            double weightSum = 0.0;
            for (unsigned int i = 0; i < nbODFs; ++i)
            {
                weightSum += weightValues[i];
                std::fill(workValue.begin(), workValue.end(), 0.0);
                anima::sphere_log_map(coefs[i], mean, workValue);
                for (unsigned int j = 0; j < m_VectorLength; ++j)
                    tangent[j] += weightValues[i] * workValue[j];
            }

            for (unsigned int j = 0; j < m_VectorLength; ++j)
                tangent[j] /= weightSum; //enables to make averaging work even if the sum of the weights is not equal to 1 (but to check)

            normTan = anima::ComputeNorm(tangent);

            std::fill(nextMean.begin(), nextMean.end(), 0.0);
            anima::sphere_exp_map(tangent, mean, nextMean);

            nIter++;
        }

        return nextMean;
    }

} //end namespace anima
