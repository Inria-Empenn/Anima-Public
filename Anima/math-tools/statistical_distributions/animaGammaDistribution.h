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
			m_ShapeParameter = 1.0;
			m_ScaleParameter = 1.0;
		}

		bool BelongsToSupport(const SingleValueType &x);
		double GetDensity(const SingleValueType &x);
		double GetLogDensity(const SingleValueType &x);
		void Fit(const MultipleValueType &sample, const std::string &method);
		void Random(MultipleValueType &sample, GeneratorType &generator);

		void SetShapeParameter(const double val);
		double GetShapeParameter() {return m_ShapeParameter;}

		void SetScaleParameter(const double val);
		double GetScaleParameter() {return m_ScaleParameter;}

	private:
		double m_ShapeParameter;
		double m_ScaleParameter;
	};
    
} // end of namespace
