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

		virtual bool BelongsToSupport(const SingleValueType &x) = 0;
		virtual double GetDensity(const SingleValueType &x) = 0;
		virtual double GetLogDensity(const SingleValueType &x) = 0;
		virtual void Fit(const MultipleValueType &sample, const std::string &method) = 0;
		virtual void Random(MultipleValueType &sample, GeneratorType &generator) = 0;
	};
} // end of namespace
