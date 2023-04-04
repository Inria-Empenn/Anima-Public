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
		virtual void Fit(const std::vector<double> &sample, const std::string &method) = 0;
		virtual void Random(std::vector<double> &sample, std::mt19937 &generator) = 0;

	private:
		double m_ShapeParameter;
		double m_ScaleParameter;
	};
} // end of namespace
