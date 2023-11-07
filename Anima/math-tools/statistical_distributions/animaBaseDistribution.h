#pragma once

#include <limits>
#include <random>
#include <string>
#include <vector>

#include <AnimaStatisticalDistributionsExport.h>

namespace anima
{
	template <typename TValueType>
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT BaseDistribution
	{
	public:
		using ValueType = TValueType;
		using SampleType = std::vector<ValueType>;
		using GeneratorType = std::mt19937;

		BaseDistribution() {}

		virtual bool BelongsToSupport(const ValueType &x) = 0;
		virtual double GetDensity(const ValueType &x) = 0;
		virtual double GetLogDensity(const ValueType &x) = 0;
		virtual void Fit(const SampleType &sample, const std::string &method) = 0;
		virtual void Random(SampleType &sample, GeneratorType &generator) = 0;
		virtual ValueType GetMean() = 0;
		virtual double GetVariance() = 0;

		double GetEpsilon() { return std::sqrt(std::numeric_limits<double>::epsilon()); }
	};
} // end of namespace
