#pragma once

#include <animaBaseDistribution.h>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT GammaDistribution : public BaseDistribution
	{
	public:
		using Distribution = std::gamma_distribution<>;
		using Generator = std::mt19937;
		double GetDensity(const double &x);
		void Fit(const std::vector<double> &sample, const std::string &method);
		void Random(std::vector<double> &sample, Generator &generator);
	};
    
} // end of namespace
