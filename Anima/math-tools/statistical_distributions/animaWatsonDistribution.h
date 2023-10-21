#pragma once

#include <animaBaseDistribution.h>

#include <itkVector.h>

namespace anima
{   
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT WatsonDistribution : public BaseDistribution<itk::Vector<double,3>,std::vector<itk::Vector<double,3>>>
	{
	public:
		using DistributionType = std::uniform_real_distribution<double>;

		WatsonDistribution()
		{
			m_MeanAxis[0] = 0;
            m_MeanAxis[1] = 0;
            m_MeanAxis[2] = 1;
            m_ConcentrationParameter = 1.0;
		}

        double GetDensity(const SingleValueType &val);
		double GetLogDensity(const SingleValueType &val);
		void Fit(const MultipleValueType &sample, const std::string &method);
		void Random(MultipleValueType &sample, GeneratorType &generator);

		void SetMeanAxis(const itk::Vector<double,3> &x);
        SingleValueType GetMeanAxis() {return m_MeanAxis;}

        void SetConcentrationParameter(const double &x);
        double GetConcentrationParameter() {return m_ConcentrationParameter;}

        void GetStandardWatsonSHCoefficients(
            std::vector<double> &coefficients, 
            std::vector<double> &derivatives
        );

    private:
        itk::Vector<double,3> m_MeanAxis;
        double m_ConcentrationParameter;
	};
    
} // end of namespace anima
