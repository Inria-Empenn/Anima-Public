#pragma once

#include <animaBaseDistribution.h>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT GammaDistribution : public BaseDistribution<double>
	{
	public:
		using DistributionType = std::gamma_distribution<double>;

		GammaDistribution()
		{
			m_ShapeParameter = 1.0;
			m_ScaleParameter = 1.0;
		}

		bool BelongsToSupport(const ValueType &x);
		double GetDensity(const ValueType &x);
		double GetLogDensity(const ValueType &x);
		void Fit(const SampleType &sample, const std::string &method);
		void Random(SampleType &sample, GeneratorType &generator);
		ValueType GetMean() { return m_ShapeParameter * m_ScaleParameter; }
		double GetVariance() { return m_ShapeParameter * m_ScaleParameter * m_ScaleParameter; }
		double GetDistance(Self *otherDistribution);

		void SetShapeParameter(const double val);
		double GetShapeParameter() { return m_ShapeParameter; }

		void SetScaleParameter(const double val);
		double GetScaleParameter() { return m_ScaleParameter; }

	private:
		double m_ShapeParameter;
		double m_ScaleParameter;
	};

} // end of namespace
