#pragma once

#include <animaBaseDistribution.h>

#include <vector>
#include <vnl/vnl_matrix.h>
#include <boost/math/distributions/beta.hpp>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT DirichletDistribution : public BaseDistribution<std::vector<double>,vnl_matrix<double>>
	{
	public:
		using UniformDistributionType = std::uniform_real_distribution<double>;
		using BetaDistributionType = boost::math::beta_distribution<double>;

		DirichletDistribution()
		{
			m_ConcentrationParameters.clear();
		}

		bool BelongsToSupport(const SingleValueType &x);
		double GetDensity(const SingleValueType &x);
		double GetLogDensity(const SingleValueType &x);
		void Fit(const MultipleValueType &sample, const std::string &method);
		void Random(MultipleValueType &sample, GeneratorType &generator);

		void SetConcentrationParameters(const std::vector<double> &val);
		std::vector<double> GetConcentrationParameters() {return m_ConcentrationParameters;}

	private:
		std::vector<double> m_ConcentrationParameters;
	};
    
} // end of namespace
