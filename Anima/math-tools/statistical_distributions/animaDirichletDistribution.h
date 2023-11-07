#pragma once

#include <animaBaseDistribution.h>

#include <vector>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_trace.h>
#include <boost/math/distributions/beta.hpp>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT DirichletDistribution : public BaseDistribution<std::vector<double>>
	{
	public:
		using UniformDistributionType = std::uniform_real_distribution<double>;
		using BetaDistributionType = boost::math::beta_distribution<double>;

		DirichletDistribution()
		{
			m_ConcentrationParameters.clear();
			m_MeanValues.clear();
			m_TotalConcentration = 0.0;
		}

		bool BelongsToSupport(const ValueType &x);
		double GetDensity(const ValueType &x);
		double GetLogDensity(const ValueType &x);
		void Fit(const SampleType &sample, const std::string &method);
		void Random(SampleType &sample, GeneratorType &generator);
		ValueType GetMean() {return m_MeanValues;}
		double GetVariance() {return vnl_trace(this->GetCovarianceMatrix());}

		void SetConcentrationParameters(const std::vector<double> &val);
		std::vector<double> GetConcentrationParameters() {return m_ConcentrationParameters;}

		vnl_matrix<double> GetCovarianceMatrix();

	private:
		std::vector<double> m_ConcentrationParameters;
		ValueType m_MeanValues;
		double m_TotalConcentration;
	};
    
} // end of namespace
