#include "animaEPGProfileIntegrands.h"

namespace anima
{

std::vector <double> EPGMonoT2Integrand::operator() (double const t)
{
    // Find value in slice profile
    using PulseIteratorType = std::vector < std::pair <double, double> >::iterator;
    double pulseProfileValue = 1.0;

    if (t <= m_SlicePulseProfile[0].first)
        pulseProfileValue = m_SlicePulseProfile[0].second;
    else if (t >= m_SlicePulseProfile.back().first)
        pulseProfileValue = m_SlicePulseProfile.back().second;
    else
    {
        PulseIteratorType it = std::lower_bound(m_SlicePulseProfile.begin(), m_SlicePulseProfile.end(), std::make_pair(t, 0.0));
        PulseIteratorType itUp = it;
        ++itUp;
        if (itUp == m_SlicePulseProfile.end())
            pulseProfileValue = it->second;
        else
        {
            double weight = (t - it->first) / (itUp->first - it->first);
            pulseProfileValue = weight * itUp->second + (1.0 - weight) * it->second;
        }
    }

    double excitationProfileValue = 1.0;

    if (t <= m_SliceExcitationProfile[0].first)
        excitationProfileValue = m_SliceExcitationProfile[0].second;
    else if (t >= m_SliceExcitationProfile.back().first)
        excitationProfileValue = m_SliceExcitationProfile.back().second;
    else
    {
        PulseIteratorType it = std::lower_bound(m_SliceExcitationProfile.begin(), m_SliceExcitationProfile.end(), std::make_pair(t, 0.0));
        PulseIteratorType itUp = it;
        ++itUp;
        if (itUp == m_SliceExcitationProfile.end())
            excitationProfileValue = it->second;
        else
        {
            double weight = (t - it->first) / (itUp->first - it->first);
            excitationProfileValue = weight * itUp->second + (1.0 - weight) * it->second;
        }
    }

    double refExcitationValue = m_SignalSimulator.GetExcitationFlipAngle();
    m_SignalSimulator.SetExcitationFlipAngle(excitationProfileValue * refExcitationValue);

    std::vector <double> epgValue = m_SignalSimulator.GetValue(m_T1Value, m_T2Value, pulseProfileValue * m_FlipAngle, 1.0);

    m_SignalSimulator.SetExcitationFlipAngle(refExcitationValue);

    return epgValue;
}

} // end namespace anima
