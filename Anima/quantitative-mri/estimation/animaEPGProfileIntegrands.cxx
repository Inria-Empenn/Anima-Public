#include "animaEPGProfileIntegrands.h"

namespace anima
{

std::vector <double> EPGMonoT2Integrand::operator() (double const t)
{
    // Find value in slice profile
    using PulseIteratorType = std::vector < std::pair <double, double> >::iterator;
    double profileValue = 1.0;

    if (t <= m_SlicePulseProfile[0].first)
        profileValue = m_SlicePulseProfile[0].second;
    else if (t >= m_SlicePulseProfile.back().first)
        profileValue = m_SlicePulseProfile.back().second;
    else
    {
        PulseIteratorType it = std::lower_bound(m_SlicePulseProfile.begin(), m_SlicePulseProfile.end(), std::make_pair(t, 0.0));
        PulseIteratorType itUp = it;
        ++itUp;
        if (itUp == m_SlicePulseProfile.end())
            profileValue = it->second;
        else
        {
            double weight = (t - it->first) / (itUp->first - it->first);
            profileValue = weight * itUp->second + (1.0 - weight) * it->second;
        }
    }

    std::vector <double> epgValue = m_SignalSimulator.GetValue(m_T1Value, m_T2Value, profileValue * m_FlipAngle, 1.0);
    return epgValue;
}

} // end namespace anima
