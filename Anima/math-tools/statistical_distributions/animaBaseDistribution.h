#pragma once

#include <string>
#include <random>

#include <AnimaStatisticalDistributionsExport.h>

namespace anima
{
	template <typename TSingleValueType, typename TMultipleValueType>
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT BaseDistribution
	{
	public:
		using SingleValueType = TSingleValueType;
		using MultipleValueType = TMultipleValueType;
		using GeneratorType = std::mt19937;

		BaseDistribution()
		{
			m_ShapeParameter = 1.0;
			m_ScaleParameter = 1.0;
		}

		void SetShapeParameter(const SingleValueType val);
		SingleValueType GetShapeParameter() {return m_ShapeParameter;}

		void SetScaleParameter(const SingleValueType val);
		double GetScaleParameter() {return m_ScaleParameter;}

		void SetConcentrationParameters(const SingleValueType val);
		SingleValueType GetConcentrationParameters() {return m_ConcentrationParameters;}

		virtual double GetDensity(const SingleValueType &x) = 0;
		virtual double GetLogDensity(const SingleValueType &x) = 0;
		virtual void Fit(const MultipleValueType &sample, const std::string &method) = 0;
		virtual void Random(MultipleValueType &sample, GeneratorType &generator) = 0;

	private:
		double m_ShapeParameter;
		double m_ScaleParameter;
		SingleValueType m_ConcentrationParameters;
	};
} // end of namespace

#include <animaBaseDistribution.hxx>