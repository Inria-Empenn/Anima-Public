#pragma once
#include <cmath>

#include "animaWatsonDistribution.h"
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>

#include <itkExceptionObject.h>
#include <itkObjectFactory.h>

#include <boost/math/special_functions/bessel.hpp>

namespace anima
{

template <class VectorType, class ScalarType>
double EvaluateWatsonPDF(const VectorType &v, const VectorType &meanAxis, const ScalarType &kappa)
{
    /************************************************************************************************
     * \fn template <class VectorType, class ScalarType>
     * 	   double
     *     EvaluateWatsonPDF(const VectorType &v,
     * 						 const VectorType &meanAxis,
     * 						 const ScalarType &kappa)
     *
     * \brief	Evaluate the Watson probability density function using the definition of
     * 			Fisher et al., Statistical Analysis of Spherical Data, 1993, p.89.
     *
     * \author	Aymeric Stamm
     * \date	July 2014
     *
     * \param	v				Sample on which evaluating the PDF.
     * \param	meanAxis        Mean axis of the Watson distribution.
     * \param	kappa		   	Concentration parameter of the Watson distribution.
     **************************************************************************************************/

    if (std::abs(kappa) < 1.0e-6)
        return 1.0 / (4.0 * M_PI);
    else if (kappa > 0)
    {
        double kappaSqrt = std::sqrt(kappa);
        double c = anima::ComputeScalarProduct(v, meanAxis);
        double inExp = kappa * (c * c - 1.0);
        return kappaSqrt * std::exp(inExp) / (4.0 * M_PI * anima::EvaluateDawsonFunctionNR(kappaSqrt));
    }
    else
    {
        double Ck = std::sqrt(-kappa / M_PI) / (2.0 * M_PI * std::erf(std::sqrt(-kappa)));
        double c = anima::ComputeScalarProduct(v, meanAxis);
        double inExp = kappa * c * c;
        return Ck * std::exp(inExp);
    }
}

template <class ScalarType>
double EvaluateWatsonPDF(const vnl_vector_fixed <ScalarType,3> &v, const vnl_vector_fixed <ScalarType,3> &meanAxis, const ScalarType &kappa)
{
    if (std::abs(v.squared_magnitude() - 1.0) > 1.0e-6 || std::abs(meanAxis.squared_magnitude() - 1.0) > 1.0e-6)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Watson distribution is on the 2-sphere.",ITK_LOCATION);

    return EvaluateWatsonPDF<vnl_vector_fixed <ScalarType,3>, ScalarType>(v,meanAxis,kappa);
}

template <class ScalarType>
double EvaluateWatsonPDF(const itk::Point <ScalarType,3> &v, const itk::Point <ScalarType,3> &meanAxis, const ScalarType &kappa)
{
    if (std::abs(v.GetVnlVector().squared_magnitude() - 1.0) > 1.0e-6 || std::abs(meanAxis.GetVnlVector().squared_magnitude() - 1.0) > 1.0e-6)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Watson distribution is on the 2-sphere.",ITK_LOCATION);

    return EvaluateWatsonPDF<itk::Point <ScalarType,3>, ScalarType>(v,meanAxis,kappa);
}

template <class ScalarType>
double EvaluateWatsonPDF(const itk::Vector <ScalarType,3> &v, const itk::Vector <ScalarType,3> &meanAxis, const ScalarType &kappa)
{
    if (std::abs(v.GetSquaredNorm() - 1.0) > 1.0e-6 || std::abs(meanAxis.GetSquaredNorm() - 1.0) > 1.0e-6)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Watson distribution is on the 2-sphere.",ITK_LOCATION);

    return EvaluateWatsonPDF<itk::Vector <ScalarType,3>, ScalarType>(v,meanAxis,kappa);
}

template <class ScalarType>
    void GetStandardWatsonSHCoefficients(const ScalarType k, std::vector<ScalarType> &coefficients, std::vector<ScalarType> &derivatives)
{
    // Computes the first 7 non-zero SH coefficients of the standard Watson PDF (multiplied by 4 M_PI).
    const unsigned int nbCoefs = 7;
    coefficients.resize(nbCoefs);
    derivatives.resize(nbCoefs);
    
    double sqrtPi = std::sqrt(M_PI);
    double kSqrt = std::sqrt(k);
    double dawsonValue = anima::EvaluateDawsonIntegral(kSqrt);
    double k2 = k * k;
    double k3 = k2 * k;
    double k4 = k3 * k;
    double k5 = k4 * k;
    double k6 = k5 * k;
    double k7 = k6 * k;
    double k8 = k7 * k;
    double k9 = k8 * k;
    double k10 = k9 * k;
    double k11 = k10 * k;
    double k12 = k11 * k;
    double k13 = k12 * k;
    
    coefficients[0] = 2.0 * sqrtPi;
    derivatives[0] = 0.0;
    
    if (k < 1.0e-4) // Handles small k values by Taylor series expansion
    {
        coefficients[1] = 8.0 * sqrtPi * k2 / 105.0;
        coefficients[2] = 32.0 * sqrtPi / std::sqrt(17.0) * k4 / 19305.0;
        coefficients[3] = 128.0 * sqrtPi * k6 / 152108775.0;
        coefficients[4] = 512.0 * sqrtPi / std::sqrt(33.0) * k8 / 94670161425.0;
        coefficients[5] = 2048.0 * sqrtPi / std::sqrt(41.0) * k10 / 488493636505875.0;
        coefficients[6] = 8192.0 * sqrtPi * k12 / 26398089922508679375.0;
        derivatives[1] = 16.0 * sqrtPi * k / 105.0;
        derivatives[2] = 128.0 * sqrtPi / std::sqrt(17.0) * k3 / 19305.0;
        derivatives[3] = 256.0 * sqrtPi * k5 / 50702925.0;
        derivatives[4] = 4096.0 * sqrtPi / std::sqrt(33.0) * k7 / 94670161425.0;
        derivatives[5] = 4096.0 * sqrtPi / std::sqrt(41.0) * k9 / 97698727301175.0;
        derivatives[6] = 32768.0 * sqrtPi * k11 / 8799363307502893125.0;
        return;
    }
    
    if (k <= 1.0)
    {
        coefficients[1] = (15.0 * kSqrt * (-21.0 + 2.0 * k) + 9.0 * (35.0 + 20.0 * k + 4.0 * k2) * dawsonValue) / (16.0 * k2);
        coefficients[2] = std::sqrt(17.0) * (3.0 * kSqrt * (-225225.0 + 30030.0 * k - 7700.0 * k2 + 248.0 * k3) + 35.0 * (19305.0 + 10296.0 * k + 2376.0 * k2 + 288.0 * k3 + 16.0 * k4) * dawsonValue) / (1024.0 * k4);
        coefficients[3] = (65.0 * kSqrt * (-540571185.0 + 78343650.0 * k - 23279256.0 * k2 + 1319472.0 * k3 - 119504.0 * k4 + 1952.0 * k5) + 1155.0 * (30421755.0 + 15872220.0 * k + 3779100.0 * k2 + 530400.0 * k3 + 46800.0 * k4 + 2496.0 * k5 + 64.0 * k6) * dawsonValue) / (32768.0 * k6);
        coefficients[4] = std::sqrt(33.0) * (17.0 * kSqrt * (-35835440515875.0 + 5394582443250.0 * k - 1690125336900.0 * k2 + 113046219000.0 * k3 - 13273088400.0 * k4 + 430017120.0 * k5 - 19972800.0 * k6 + 198272.0 * k7) + 6435.0 * (94670161425.0 + 48862018800.0 * k + 11794280400.0 * k2 + 1747300800.0 * k3 + 174730080.0 * k4 + 12155136.0 * k5 + 578816.0 * k6 + 17408.0 * k7 + 256.0 * k8) * dawsonValue) / (4194304.0 * k8);
        coefficients[5] = std::sqrt(41.0) * (15.0 * kSqrt * (-1504202171771324025.0 + 231415718734049850.0 * k - 74497792937808240.0 * k2 + 5400692878220640.0 * k3 - 702763957415520.0 * k4 + 28363135120320.0 * k5 - 1839815028480.0 * k6 + 39036839424.0 * k7 - 1101379840.0 * k8 + 7371264.0 * k9) + 46189.0 * (488493636505875.0 + 250509557182500.0 * k + 60934757152500.0 * k2 + 9285296328000.0 * k3 + 984804156000.0 * k4 + 76242902400.0 * k5 + 4381776000.0 * k6 + 185472000.0 * k7 + 5564160.0 * k8 + 107520.0 * k9 + 1024.0 * k10) * dawsonValue) / (134217728.0 * k10);
        coefficients[6] = 7.0 * (15.0 * kSqrt * (-169963222029741381866625.0 + 26519084288328442560750.0 * k - 8678973039816581201700.0 * k2 + 659973586637953387800.0 * k3 - 90753965143384102560.0 * k4 + 4095016957109856960.0 * k5 - 305795964991870080.0 * k6 + 8403995054058240.0 * k7 - 349108244808960.0 * k8 + 5257552774656.0 * k9 - 99504444416.0 * k10 + 480360448.0 * k11) + 676039.0 * (3771155703215525625.0 + 1925696529301545000.0 * k + 470725818273711000.0 * k2 + 72980747019180000.0 * k3 + 8010081989910000.0 * k4 + 657237496608000.0 * k5 + 41447409696000.0 * k6 + 2030077209600.0 * k7 + 76896864000.0 * k8 + 2204928000.0 * k9 + 45619200.0 * k10 + 614400.0 * k11 + 4096.0 * k12) * dawsonValue) / (8589934592.0 * k12);
    }
    else
    {
        coefficients[1] = (15.0 * kSqrt * (-21.0 / k2 + 2.0 / k) + 9.0 * (35.0 / k2 + 20.0 / k + 4.0) * dawsonValue) / 16.0;
        coefficients[2] = std::sqrt(17.0) * (3.0 * kSqrt * (-225225.0 / k4 + 30030.0 / k3 - 7700.0 / k2 + 248.0 / k) + 35.0 * (19305.0 / k4 + 10296.0 / k3 + 2376.0 / k2 + 288.0 / k + 16.0) * dawsonValue) / 1024.0;
        coefficients[3] = (65.0 * kSqrt * (-540571185.0 / k6 + 78343650.0 / k5 - 23279256.0 / k4 + 1319472.0 / k3 - 119504.0 / k2 + 1952.0 / k) + 1155.0 * (30421755.0 / k6 + 15872220.0 / k5 + 3779100.0 / k4 + 530400.0 / k3 + 46800.0 / k2 + 2496.0 / k + 64.0) * dawsonValue) / 32768.0;
        coefficients[4] = std::sqrt(33.0) * (17.0 * kSqrt * (-35835440515875.0 / k8 + 5394582443250.0 / k7 - 1690125336900.0 / k6 + 113046219000.0 / k5 - 13273088400.0 / k4 + 430017120.0 / k3 - 19972800.0 / k2 + 198272.0 / k) + 6435.0 * (94670161425.0 / k8 + 48862018800.0 / k7 + 11794280400.0 / k6 + 1747300800.0 / k5 + 174730080.0 / k4 + 12155136.0 / k3 + 578816.0 / k2 + 17408.0 / k + 256.0) * dawsonValue) / 4194304.0;
        coefficients[5] = std::sqrt(41.0) * (15.0 * kSqrt * (-1504202171771324025.0 / k10 + 231415718734049850.0 / k9 - 74497792937808240.0 / k8 + 5400692878220640.0 / k7 - 702763957415520.0 /  k6 + 28363135120320.0 / k5 - 1839815028480.0 / k4 + 39036839424.0 / k3 - 1101379840.0 / k2 + 7371264.0 / k) + 46189.0 * (488493636505875.0 / k10 + 250509557182500.0 / k9 + 60934757152500.0 / k8 + 9285296328000.0 / k7 + 984804156000.0 / k6 + 76242902400.0 / k5 + 4381776000.0 / k4 + 185472000.0 / k3 + 5564160.0 / k2 + 107520.0 / k + 1024.0) * dawsonValue) / 134217728.0;
        coefficients[6] = 7.0 * (15.0 * kSqrt * (-169963222029741381866625.0 / k12 + 26519084288328442560750.0 / k11 - 8678973039816581201700.0 / k10 + 659973586637953387800.0 / k9 - 90753965143384102560.0 / k8 + 4095016957109856960.0 / k7 - 305795964991870080.0 / k6 + 8403995054058240.0 / k5 - 349108244808960.0 / k4 + 5257552774656.0 / k3 - 99504444416.0 / k2 + 480360448.0 / k) + 676039.0 * (3771155703215525625.0 / k12 + 1925696529301545000.0 / k11 + 470725818273711000.0 / k10 + 72980747019180000.0 / k9 + 8010081989910000.0 / k8 + 657237496608000.0 / k7 + 41447409696000.0 / k6 + 2030077209600.0 / k5 + 76896864000.0 / k4 + 2204928000.0 / k3 + 45619200.0 / k2 + 614400.0 / k + 4096.0) * dawsonValue) / 8589934592.0;
    }
    
    for (unsigned int i = 1;i < nbCoefs;++i)
        coefficients[i] *= (sqrtPi / dawsonValue);
    
    // checked
    double tmpValue1 = -21.0 + 2.0 * k;
    double tmpValue2 = 35.0 + 20.0 * k + 4.0 * k2;
    derivatives[1] = sqrtPi * (60.0 * kSqrt + 15.0 * tmpValue1 / kSqrt + 30.0 * kSqrt * tmpValue1 + 9.0 * tmpValue2 / kSqrt + 9.0 * (20.0 + 8.0 * k) * 2.0 * dawsonValue) / (16.0 * k2 * 2.0 * dawsonValue) - sqrtPi * (30.0 * kSqrt * tmpValue1 + 9.0 * tmpValue2 * 2.0 * dawsonValue) / (16.0 * k2 * kSqrt * 4.0 * dawsonValue * dawsonValue) - sqrtPi * (30.0 * kSqrt * tmpValue1 + 9.0 * tmpValue2 * 2.0 * dawsonValue) / (8.0 * k3 * 2.0 * dawsonValue);
    
    // checked
    tmpValue1 = -225225.0 + 30030.0 * k - 7700.0 * k2 + 248.0 * k3;
    tmpValue2 = 19305.0 + 10296.0 * k + 2376.0 * k2 + 288.0 * k3 + 16.0 * k4;
    derivatives[2] = std::sqrt(17.0) * sqrtPi * (6.0 * kSqrt * (30030.0 - 15400.0 * k + 744.0 * k2) + 3.0 * tmpValue1 / kSqrt + 6.0 * kSqrt * tmpValue1 + 35.0 * tmpValue2 / kSqrt + 35.0 * (10296.0 + 4752.0 * k + 864.0 * k2 + 64.0 * k3) * 2.0 * dawsonValue) / (1024.0 * k4 * 2.0 * dawsonValue) - std::sqrt(17.0) * sqrtPi * (6.0 * kSqrt * tmpValue1 + 35.0 * tmpValue2 * 2.0 * dawsonValue) / (1024.0 * k4 * kSqrt * 4.0 * dawsonValue * dawsonValue) - std::sqrt(17.0) * sqrtPi * (6.0 * kSqrt * tmpValue1 + 35.0 * tmpValue2 * 2.0 * dawsonValue) / (256.0 * k5 * 2.0 * dawsonValue);
    
    // checked
    tmpValue1 = -540571185.0 + 78343650.0 * k - 23279256.0 * k2 + 1319472.0 * k3 - 119504.0 * k4 + 1952.0 * k5;
    tmpValue2 = 30421755.0 + 15872220.0 * k + 3779100.0 * k2 + 530400.0 * k3 + 46800.0 * k4 + 2496.0 * k5 + 64.0 * k6;
    derivatives[3] = sqrtPi * (130.0 * kSqrt * (78343650.0 - 46558512.0 * k + 3958416.0 * k2 - 478016.0 * k3 + 9760.0 * k4) + 65.0 * tmpValue1 / kSqrt + 130.0 * kSqrt * tmpValue1 + 1155.0 * tmpValue2 / kSqrt + 1155.0 * (15872220.0 + 7558200.0 * k + 1591200.0 * k2 + 187200.0 * k3 + 12480.0 * k4 + 384.0 * k5) * 2.0 * dawsonValue) / (32768.0 * k6 * 2.0 * dawsonValue) - sqrtPi * (130.0 * kSqrt * tmpValue1 + 1155.0 * tmpValue2 * 2.0 * dawsonValue) / (32768.0 * k6 * kSqrt * 4.0 * dawsonValue * dawsonValue) - 3.0 * sqrtPi * (130.0 * kSqrt * tmpValue1 + 1155.0 * tmpValue2 * 2.0 * dawsonValue) / (16384.0 * k7 * 2.0 * dawsonValue);
    
    // checked
    tmpValue1 = -35835440515875.0 + 5394582443250.0 * k - 1690125336900.0 * k2 + 113046219000.0 * k3 - 13273088400.0 * k4 + 430017120.0 * k5 - 19972800.0 * k6 + 198272.0 * k7;
    tmpValue2 = 94670161425.0 + 48862018800.0 * k + 11794280400.0 * k2 + 1747300800.0 * k3 + 174730080.0 * k4 + 12155136.0 * k5 + 578816.0 * k6 + 17408.0 * k7 + 256.0 * k8;
    derivatives[4] = std::sqrt(33.0) * sqrtPi * (34.0 * kSqrt * (5394582443250.0 - 3380250673800.0 * k + 339138657000.0 * k2 - 53092353600.0 * k3 + 2150085600.0 * k4 - 119836800.0 * k5 + 1387904.0 * k6) + 17.0 * tmpValue1 / kSqrt + 34.0 * kSqrt * tmpValue1 + 6435.0 * tmpValue2 / kSqrt + 6435.0 * (48862018800.0 + 23588560800.0 * k + 5241902400.0 * k2 + 698920320.0 * k3 + 60775680.0 * k4 + 3472896.0 * k5 + 121856.0 * k6 + 2048.0 * k7) * 2.0 * dawsonValue) / (4194304.0 * k8 * 2.0 * dawsonValue) - std::sqrt(33.0) * sqrtPi * (34.0 * kSqrt * tmpValue1 + 6435.0 * tmpValue2 * 2.0 * dawsonValue) / (4194304.0 * k8 * kSqrt * 4.0 * dawsonValue * dawsonValue) - std::sqrt(33.0) * sqrtPi * (34.0 * kSqrt * tmpValue1 + 6435.0 * tmpValue2 * 2.0 * dawsonValue) / (524288.0 * k9 * 2.0 * dawsonValue);
    
    // checked
    tmpValue1 = -1504202171771324025.0 + 231415718734049850.0 * k - 74497792937808240.0 * k2 + 5400692878220640.0 * k3 - 702763957415520.0 * k4 + 28363135120320.0 * k5 - 1839815028480.0 * k6 + 39036839424.0 * k7 - 1101379840.0 * k8 + 7371264.0 * k9;
    tmpValue2 = 488493636505875.0 + 250509557182500.0 * k + 60934757152500.0 * k2 + 9285296328000.0 * k3 + 984804156000.0 * k4 + 76242902400.0 * k5 + 4381776000.0 * k6 + 185472000.0 * k7 + 5564160.0 * k8 + 107520.0 * k9 + 1024.0 * k10;
    derivatives[5] = std::sqrt(41.0) * sqrtPi * (30.0 * kSqrt * (231415718734049850.0 - 148995585875616480.0 * k + 16202078634661920.0 * k2 - 2811055829662080.0 * k3 + 141815675601600.0 * k4 - 11038890170880.0 * k5 + 273257875968.0 * k6 - 8811038720.0 * k7 + 66341376.0 * k8) + 15.0 * tmpValue1 / kSqrt + 30.0 * kSqrt * tmpValue1 + 46189.0 * tmpValue2 / kSqrt + 46189.0 * (250509557182500.0 + 121869514305000.0 * k + 27855888984000.0 * k2 + 3939216624000.0 * k3 + 381214512000.0 * k4 + 26290656000.0 * k5 + 1298304000.0 * k6 + 44513280.0 * k7 + 967680.0 * k8 + 10240.0 * k9) * 2.0 * dawsonValue) / (134217728.0 * k10 * 2.0 * dawsonValue) - std::sqrt(41.0) * sqrtPi * (30.0 * kSqrt * tmpValue1 + 46189.0 * tmpValue2 * 2.0 * dawsonValue) / (134217728.0 * k10 * kSqrt * 4.0 * dawsonValue * dawsonValue) - 5.0 * std::sqrt(41.0) * sqrtPi * (30.0 * kSqrt * tmpValue1 + 46189.0 * tmpValue2 * 2.0 * dawsonValue) / (67108864.0 * k11 * 2.0 * dawsonValue);
    
    // checked
    tmpValue1 = -169963222029741381866625.0 + 26519084288328442560750.0 * k - 8678973039816581201700.0 * k2 + 659973586637953387800.0 * k3 - 90753965143384102560.0 * k4 + 4095016957109856960.0 * k5 - 305795964991870080.0 * k6 + 8403995054058240.0 * k7 - 349108244808960.0 * k8 + 5257552774656.0 * k9 - 99504444416.0 * k10 + 480360448.0 * k11;
    tmpValue2 = 3771155703215525625.0 + 1925696529301545000.0 * k + 470725818273711000.0 * k2 + 72980747019180000.0 * k3 + 8010081989910000.0 * k4 + 657237496608000.0 * k5 + 41447409696000.0 * k6 + 2030077209600.0 * k7 + 76896864000.0 * k8 + 2204928000.0 * k9 + 45619200.0 * k10 + 614400.0 * k11 + 4096.0 * k12;
    derivatives[6] = 7.0 * sqrtPi * (30.0 * kSqrt * (26519084288328442560750.0 - 17357946079633162403400.0 * k + 1979920759913860163400.0 * k2 - 363015860573536410240.0 * k3 + 20475084785549284800.0 * k4 - 1834775789951220480.0 * k5 + 58827965378407680.0 * k6 - 2792865958471680.0 * k7 + 47317974971904.0 * k8 - 995044444160.0 * k9 + 5283964928.0 * k10) + 15.0 * tmpValue1 / kSqrt + 30.0 * kSqrt * tmpValue1 + 676039.0 * tmpValue2 / kSqrt + 676039.0 * (1925696529301545000.0 + 941451636547422000.0 * k + 218942241057540000.0 * k2 + 32040327959640000.0 * k3 + 3286187483040000.0 * k4 + 248684458176000.0 * k5 + 14210540467200.0 * k6 + 615174912000.0 * k7 + 19844352000.0 * k8 + 456192000.0 * k9 + 6758400.0 * k10 + 49152.0 * k11) * 2.0 * dawsonValue) / (8589934592.0 * k12 * 2.0 * dawsonValue) - 7.0 * sqrtPi * (30.0 * kSqrt * tmpValue1 + 676039.0 * tmpValue2 * 2.0 * dawsonValue) / (8589934592.0 * k12 * kSqrt * 4.0 * dawsonValue * dawsonValue) - 21.0 * sqrtPi * (30.0 * kSqrt * tmpValue1 + 676039.0 * tmpValue2 * 2.0 * dawsonValue) / (2147483648.0 * k13 * 2.0 * dawsonValue);
}

} // end namespace anima
