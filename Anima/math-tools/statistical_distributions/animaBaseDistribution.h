#pragma once

#include <string>
#include <vector>
#include <random>

#include <AnimaStatisticalDistributionsExport.h>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT BaseDistribution
	{
	public:
		using GeneratorType = std::mt19937;
		using VectorType = std::vector<double>;

		BaseDistribution()
		{
			m_ShapeParameter = 1.0;
			m_ScaleParameter = 1.0;
		}

		void SetShapeParameter(const double val) {m_ShapeParameter = val;}
		double GetShapeParameter() {return m_ShapeParameter;}

		void SetScaleParameter(const double val) {m_ScaleParameter = val;}
		double GetScaleParameter() {return m_ScaleParameter;}

		virtual double GetDensity(const double &x) = 0;
		virtual double GetLogDensity(const double &x);
		virtual void Fit(const VectorType &sample, const std::string &method) = 0;
		virtual void Random(VectorType &sample, GeneratorType &generator) = 0;

	private:
		double m_ShapeParameter;
		double m_ScaleParameter;
	};
} // end of namespace
