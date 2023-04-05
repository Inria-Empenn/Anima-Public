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
		double GetDensity(const SingleValueType &x);
		double GetLogDensity(const SingleValueType &x);
		void Fit(const MultipleValueType &sample, const std::string &method);
		void Random(MultipleValueType &sample, GeneratorType &generator);
	};
    
} // end of namespace
