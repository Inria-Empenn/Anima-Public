#pragma once

#include <animaBaseDistribution.h>

#include <vector>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT GammaDistribution : public BaseDistribution<double,std::vector<double>>
	{
	public:
		using DistributionType = std::gamma_distribution<double>;

		GammaDistribution()
		{
			this->SetShapeParameter(1.0);
			this->SetScaleParameter(1.0);
		}

		void SetShapeParameter(const SingleValueType val);
		void SetScaleParameter(const SingleValueType val);
		double GetDensity(const SingleValueType &x);
		double GetLogDensity(const SingleValueType &x);
		void Fit(const MultipleValueType &sample, const std::string &method);
		void Random(MultipleValueType &sample, GeneratorType &generator);
	};
    
} // end of namespace
