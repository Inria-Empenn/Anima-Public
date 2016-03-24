#pragma once

#include "AnimaGraphCutSegmentationExport.h"
#include "animaModelInitializer.h"
#include "time.h"

namespace anima
{

/** @brief Class initializing ramdomly a gaussian model
 *
   */
class ANIMAGRAPHCUTSEGMENTATION_EXPORT RandomInitializer: public ModelInitializer
{
public:

    /** Standard class typedefs. */
    typedef RandomInitializer  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(RandomInitializer, itk::ProcessObject);

    typedef itk::VariableLengthVector<double> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;

    void Update();

    /** @brief Set min values of random means
       */
    void SetMinValues(std::vector<double> &min){this->minValues=min;}

    /** @brief Set max values of random means
       */
    void SetMaxValues(std::vector<double> &max){this->maxValues=max;}

    itkSetMacro(NbGaussian, unsigned int);
    itkGetMacro(NbGaussian, unsigned int);

    itkSetMacro(DimensionGaussian, unsigned int);
    itkGetMacro(DimensionGaussian, unsigned int);

protected:

    RandomInitializer(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    RandomInitializer()
    {
        srand(time(NULL));

        m_NbGaussian = 1;
        m_DimensionGaussian= 1;

    }
    virtual ~RandomInitializer(){}

    /** @brief initilize a Gaussian distribution
     *
       */
    itk::Statistics::GaussianMembershipFunction<itk::VariableLengthVector<double> >::Pointer randomDistribution();

    double randUniform(double min,double max);


    /** @brief min values
        */
    std::vector<double> minValues;

    /** @brief max values
       */
    std::vector<double> maxValues;

    unsigned int m_NbGaussian;
    unsigned int m_DimensionGaussian;

};

}

