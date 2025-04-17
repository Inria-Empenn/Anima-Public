#include "animaODFAverageImageFilter.h"
#include <animaLogExpMapsUnitSphere.h>
#include <animaVectorOperations.h>

#include <itkImageRegionIterator.h>


namespace anima
{
    /**
     * animaODFAverage computes the average of input ODFs voxel by voxel by using a Riemannian framework. 
     * Here we denote by ODF either a diffusion ODF, a fiber ODF or a TOD (see Dhollander et al., 2014 for the latter).
     * This code has been deviced only for situations where at most 1 input ODF is not normalized (normalized = integrates to 1 over the sphere).
     * It's the case when we compute an average between several TOD images (all ODF are normalized), or an average between a fiber ODF image and a TOD image (only the ODFs from the fiber ODF image are not normalized), which are the only two situations necessary for now.
     * But using for example this code to compute the average between 2 fiber ODF images doesn't enter in the scope of the following algorithm (even if an output image can still be computed), as then more than 1 ODF is not normalized (in all voxels).
     */

    ODFAverageImageFilter::ODFAverageImageFilter()
    {
        m_EpsValueTestNormalized = std::pow(10, -6);
        m_EpsValueTestNull = std::pow(10, -12);
        m_VectorLength = 0;
        m_ODFSHBasis = nullptr;
        m_SamplePoints = g_SphericalDesignPoints;
        //Put sample points (ie spherical designs) in spherical coordinates
        std::for_each(m_SamplePoints.begin(), m_SamplePoints.end(), [](VectorType &v){anima::TransformCartesianToSphericalCoordinates(v,v);});
    }


    ODFAverageImageFilter::~ODFAverageImageFilter() = default;


    void ODFAverageImageFilter::AddWeightImage(const unsigned int i, const WeightImagePointer &weightImage)
    {
        /**
         * Adds a newly read weight image to the vector m_WeightImages 
         */
        if (i == m_WeightImages.size())
        {
            m_WeightImages.push_back(weightImage);
            return;
        }

        if (i > m_WeightImages.size())
        {
            itkExceptionMacro("Weight images must be added contiguously.");
        }

        m_WeightImages[i] = weightImage;
    }


    void ODFAverageImageFilter::BeforeThreadedGenerateData()
    {
        /**
         * Preliminary computations
         */
        Superclass::BeforeThreadedGenerateData();

        m_VectorLength = this->GetInput(0)->GetVectorLength();
        //Get the maximum spherical harmonic (SH) order from the vector length
        unsigned int odfSHOrder = std::round(-1.5 + 0.5 * std::sqrt(8.0 * static_cast<double>(m_VectorLength) + 1.0));

        //Build m_SpherHarm (it will contain the values of each spherical harmonic (columns) for each sample point (lines))
        m_SpherHarm.set_size(NB_SAMPLE_POINTS, m_VectorLength);
        m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(odfSHOrder);
        double sqrt2 = std::sqrt(2.0);
        VectorType currentPoint;
        for (unsigned int i = 0; i < NB_SAMPLE_POINTS; i++)
        {
            currentPoint = m_SamplePoints[i];
            unsigned int c = 0;
            for (double l = 0; l <= odfSHOrder; l += 2)
                for (double m = -l; m <= l; m++)
                    m_SpherHarm.put(i, c++, m_ODFSHBasis->getNthSHValueAtPosition(l, m, currentPoint[0], currentPoint[1])); 
        }

        //Pseudo-inverse computation (will be used to move from sampled representation to representation by coefficients)
        m_SolveSHMatrix = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose();
    }


    void ODFAverageImageFilter::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
    {
        /**
         *Main function. For each voxel of the output image, it computes the average of the ODFs from the corresponding input voxels in input images. 
         */

        //Initializations
        using InputImageIteratorType = itk::ImageRegionConstIterator<InputImageType>;
        using InputWeightIteratorType = itk::ImageRegionConstIterator<WeightImageType>;
        using OutputImageIteratorType = itk::ImageRegionIterator<OutputImageType>;

        unsigned int numInputs = this->GetNumberOfIndexedInputs(); //number of input ODFs

        std::vector<InputImageIteratorType> inItrs(numInputs); //vector of input images (which are here iterators on the input images)
        std::vector<InputWeightIteratorType> weightItrs(numInputs); //vector of iterators on weight images
        for (unsigned int i = 0; i < numInputs; ++i)
        {
            inItrs[i] = InputImageIteratorType(this->GetInput(i), outputRegionForThread);
            weightItrs[i] = InputWeightIteratorType(m_WeightImages[i], outputRegionForThread);
        }
        OutputImageIteratorType outItr(this->GetOutput(), outputRegionForThread); //iterator on output image

        VectorType coefODFValue(m_VectorLength);
        VectorType sampledODFValue(NB_SAMPLE_POINTS);
        Vector2DType currentODFs(numInputs);
        Vector2DType notNullCurrentODFs;
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

            //Compute sum of weights of all input images for the current voxel
            for (unsigned int i = 0; i < numInputs; ++i)
            {
                double weightValue = weightItrs[i].Get();
                weightSum += weightValue;
                currentWeights[i] = weightValue;
            }

            if (weightSum < m_EpsValueTestNull)
            //If in a voxel, the sum of weights is 0, then we output a null distribution
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
                this->IncrementNumberOfProcessedPoints(); //to see the progression
            }

            else
            {
                for (unsigned int i = 0; i < numInputs; ++i)
                {
                    //Read input ODFs in the current voxel position and change their type 
                    this->ConvertODF(inItrs[i].Get(), inputValue);
                    currentODFs[i] = inputValue;
                }

                //Remove null ODFs and their associated weights to avoid useless computations
                numNotNullODFs = this->RemoveNullODFs(currentODFs, currentWeights, numInputs, notNullCurrentODFs, currentWeightsNotNullODFs);
                switch (numNotNullODFs)
                {
                    case 0:
                        //Here, all input ODFs are null, so we simply output a null distribution
                        outputValue.Fill(0.0);
                        break;
                    
                    case 1:
                        //Here, only 1 input ODF is not null, so we simply output this ODF
                        for (unsigned int i = 0; i < m_VectorLength; i++)
                        {
                            outputValue[i] = notNullCurrentODFs[0][i];
                        }
                        break;
                    
                    default:
                        //Here, more than 1 ODF is not null, so we need to compute an average with the Riemannian framework
                        indexODFtoNormalize = this->GetIndexODFtoNormalize(notNullCurrentODFs, numNotNullODFs);
                        for (unsigned int i = 0; i < numNotNullODFs; i++)
                        {
                            this->DiscretizeODF(notNullCurrentODFs[i], sampledODFValue); //it enables to compute easily the square root
                            coefODFValue = this->GetSquareRootODFCoef(sampledODFValue); //the averaging needs to be performed on square roots of ODFs (see Goh et al., 2011)
                            if (i == indexODFtoNormalize)
                            {
                                //Here, the current ODF needs to be normalized (it's the only one)
                                normValue = this->NormalizeODF(coefODFValue, notNullCurrentODFs[i]);
                            }
                            else
                            {
                                notNullCurrentODFs[i] = coefODFValue;
                            }
                        }
                        //notNullCurrentODFs now contains the ODFs in the correct form for averaging
                        coefODFValue = this->ComputeAverage(notNullCurrentODFs, currentWeightsNotNullODFs, numNotNullODFs);
                        if (indexODFtoNormalize != -1)
                        {
                            //Here, one input ODF has been normalized, so we reverse the normalization to retrieve the initial ODF amplitude
                            this->InverseNormalization(normValue, coefODFValue);
                        }
                        this->DiscretizeODF(coefODFValue, sampledODFValue); //it enables to compute easily the square
                        outputValue = this->GetSquareODFCoef(sampledODFValue); //reverse the square root operation
                }

                //Reset vectors
                notNullCurrentODFs.clear();
                currentWeightsNotNullODFs.clear();

                //Set output value and move to next voxel
                outItr.Set(outputValue);

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
         * Deletes pointers which are not SmartPointers
         */
        delete m_ODFSHBasis;
    }


    void ODFAverageImageFilter::ConvertODF(const InputPixelType &inputODF, VectorType &convertedODF)
    {
        /**
         * Converts an ODF from InputPixelType to VectorType
         */
        for (unsigned int i = 0; i < m_VectorLength; i++)
        {
            convertedODF[i] = inputODF[i]; 
        }
    }


    int ODFAverageImageFilter::RemoveNullODFs(const Vector2DType &allODFs, const VectorType &allWeights, int nbODFs, Vector2DType &selectedODFs, VectorType &selectedWeights)
    {
        /**
         * allODFs contains all input ODFs in a given voxel position and allWeights contains the associated weights. nbODFs is the number of input ODFs.
         * selectedODFs will contain only the input ODFs which are not null, and we will keep only their associated weights in selectedWeights.
         * The function also returns the number of not null ODFs.
        */
        int numNotNullODFs = 0;
        for (unsigned int i = 0; i < nbODFs; i++)
        {
            VectorType odf = allODFs[i];

            //Check is the ODF is null (namely all its coefficients are 0)
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

            //Keep the ODF if it's not null
            if (!currentODFisNull) 
            {
                selectedODFs.push_back(odf);
                selectedWeights.push_back(allWeights[i]);
                numNotNullODFs++;
            }
        }
        return numNotNullODFs;
    }


    int ODFAverageImageFilter::GetIndexODFtoNormalize(Vector2DType &odfs, int &numNotNullODFs)
    {
        /**
         * Gets the position of the potential ODF that needs to be normalized inside the vector of ODFs odfs.
         * Returns -1 if all input ODFs are already normalized.
         * 
         */
        int index = -1;
        int i=0;
        while (index == -1 && i < numNotNullODFs)
        {
            if (std::abs(odfs[i][0]-1/std::sqrt(4*M_PI)) >= m_EpsValueTestNormalized) 
            //An ODF is normalized if and only if its first coefficient is equal to 1/sqrt(4*PI)
            {
                index = i;
            }
            i++;
        }
        return index;

    }


    double ODFAverageImageFilter::NormalizeODF(const VectorType &inputODF, VectorType &normalizedODF)
    {
        /**
         * Normalizes an ODF so that its square integrates to 1 over the sphere
         * (This function is applied to square root of input ODFs, so the aim is now to have a square integral equal to 1)
         */
        double sumSquareCoefs = 0;
        double normValue;

        //Compute the square integral of the ODF (which is simply the square sum of its coefficients, since SH form an orthornormal basis)
        for (unsigned int i = 0; i < m_VectorLength; i++)
        {
            sumSquareCoefs+=std::pow(inputODF[i],2);
        }
        //If the square sum is null, we output null vector
        if (std::abs(sumSquareCoefs)<m_EpsValueTestNull)
        {
            normValue = 0.0;
            for (unsigned int i = 0; i < m_VectorLength; i++)
            {
                normalizedODF[i] = normValue;
            }
        }
        else
        {
            //Multiply each coefficient by normValue so that the square integral becomes 1 (ie the ODF becomes normalized)
            normValue = 1/std::sqrt(sumSquareCoefs);
            for (unsigned int i = 0; i < m_VectorLength; i++)
            {
                normalizedODF[i] = inputODF[i]*normValue;
            }
        }
        return normValue;
    }


    void ODFAverageImageFilter::DiscretizeODF(const VectorType &coefODF, VectorType &sampledODF)
    {
        /**
         * Moves an ODF from representation by SH coefficients to a sample form (namely vector with values at sample points)
         */
        VectorType currentPoint;
        for (unsigned int i = 0; i < NB_SAMPLE_POINTS; i++)
        {
            currentPoint = m_SamplePoints[i];
            sampledODF[i] = m_ODFSHBasis->getValueAtPosition(coefODF, currentPoint[0], currentPoint[1]); //currentPoint[0] : theta coordinate. currentPoint[1] : phi coordinate
        }
    }


    void ODFAverageImageFilter::InverseNormalization(const double &normValue, VectorType &odf)
    {
        /**If 1 (and so exactly 1) input ODF has been normalized by multiplication by a scalar normValue, we here divide the output average ODF by normValue to retrieve the initial ODF amplitude.
         * In our case, if an input ODF has been normalized, it's the fiber ODF before an averaging with a TOD. And indeed, in that case, the amplitudes of the output average ODF and its peaks must be the ones of the input fiber ODF.
        */
        if (!std::abs(normValue) < m_EpsValueTestNull)
        {
            for (unsigned int j = 0; j < m_VectorLength; j++)
            {
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
         */

        MatrixType squareRootOdf(NB_SAMPLE_POINTS, 1);
        MatrixType squareRootCoef(m_VectorLength, 1);

        //Compute the square root of each sample point
        for (unsigned int i = 0; i < NB_SAMPLE_POINTS; ++i)
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
         * Returns the square of an input ODF
         * Output ODF is represented by coefficients
         * Input ODF is in sample form
         */

        MatrixType squareOdf(NB_SAMPLE_POINTS, 1);
        MatrixType squareCoef(m_VectorLength, 1);

        //Compute the square of each sample point
        for (unsigned int i = 0; i < NB_SAMPLE_POINTS; ++i)
            squareOdf.put(i, 0, std::pow(odf[i], 2));

        //Convert from sample form to representation by coefficients
        squareCoef = m_SolveSHMatrix * squareOdf;

        //Store the result in a vector instead of a matrix
        OutputPixelType squareModelValue(m_VectorLength);
        for (unsigned int i = 0; i < m_VectorLength; ++i)
            squareModelValue[i] = squareCoef.get(i, 0);

        return squareModelValue;
    }


    ODFAverageImageFilter::VectorType ODFAverageImageFilter::ComputeAverage(const Vector2DType &coefs, const VectorType &weightValues, int nbODFs)
    {
        /**
         * Computes the weighted average of input ODFs using the method of Goh et al., 2O11 (Algorithm 1)
         * coefs: Vector of size nbODFs, each element being a vector of size m_VectorLength representing an ODF
         * weightValues: Vector of size nbODFs with the weights of the input ODFs
         * Returns the average ODF represented by coefficients
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
                tangent[j] /= weightSum; //enables to make averaging work even if the sum of the weights is not equal to 1
            normTan = anima::ComputeNorm(tangent);

            std::fill(nextMean.begin(), nextMean.end(), 0.0);
            anima::sphere_exp_map(tangent, mean, nextMean);

            nIter++;
        }

        return nextMean;
    }
} //end namespace anima
