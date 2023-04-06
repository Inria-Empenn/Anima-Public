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

		BaseDistribution() {}

		virtual void SetShapeParameter(const SingleValueType val) {m_ShapeParameter = val;}
		SingleValueType GetShapeParameter() {return m_ShapeParameter;}

		virtual void SetScaleParameter(const SingleValueType val) {m_ScaleParameter = val;}
		SingleValueType GetScaleParameter() {return m_ScaleParameter;}

		virtual void SetConcentrationParameter(const SingleValueType val) {m_ConcentrationParameter = val;}
		SingleValueType GetConcentrationParameter() {return m_ConcentrationParameter;}

		virtual double GetDensity(const SingleValueType &x) = 0;
		virtual double GetLogDensity(const SingleValueType &x) = 0;
		virtual void Fit(const MultipleValueType &sample, const std::string &method) = 0;
		virtual void Random(MultipleValueType &sample, GeneratorType &generator) = 0;

	private:
		SingleValueType m_ShapeParameter;
		SingleValueType m_ScaleParameter;
		SingleValueType m_ConcentrationParameter;
	};
} // end of namespace
