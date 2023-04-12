#pragma once

#include "animaDDIDistribution.h"

#include <animaVectorOperations.h>
#include <animaTrigonometricFunctions.h>
#include <animaHyperbolicFunctions.h>
#include <animaBesselFunctions.h>

#include <itkTimeProbe.h>

namespace anima
{

template <class ScalarType> double ComputeAsymmetricPDF(const std::vector<ScalarType> &w, const std::vector<ScalarType> &mu,
                                                        const ScalarType &kappa, const ScalarType &d, const ScalarType &nu,
                                                        const ScalarType &step, std::vector<ScalarType> &integrand)
{
    double invNu = 1.0 - nu;

    if (invNu < 0.0)
        invNu = 0.0;

    if (invNu > 1.0)
        invNu = 1.0;

    if (invNu < 1e-4 || d < 1e-16)
        return 0;

    double rNu = nu / invNu;
    double gaussDenom = 2.0*invNu*d;
    double C = (kappa + 1.0) / std::pow(gaussDenom*M_PI,1.5);
    double wNorm = anima::ComputeNorm(w);
    double sphereRadius = std::sqrt(nu * d);

    double resVal;

    if (nu < 1e-8)
    {
        double wpara = anima::ComputeScalarProduct(mu, w);
        double wparaSq = wpara * wpara;
        double wperpSq = wNorm * wNorm - wparaSq;
        resVal = C * exp(-((kappa + 1.0)*wperpSq+wparaSq)/gaussDenom);
    }
    else
    {
        if (kappa < 1e-4)
        {
            wNorm /= sphereRadius;

            C *= ((1.0 - std::exp(-2.0*rNu*wNorm)) / (2.0*rNu*wNorm));
            resVal = C * std::exp(-rNu*wNorm*wNorm/2.0
                                  +rNu*wNorm
                                  -rNu/2.0);
        }
        else
        {
            C *= (kappa / (1.0 - exp(-2.0*kappa)));

            double wpara = anima::ComputeScalarProduct(mu, w);
            double wparaSq = wpara * wpara;
            double wperpSq = wNorm * wNorm - wparaSq;

            if (wperpSq < 0)
                wperpSq = 0;

            double wperp = std::sqrt(wperpSq);

            // Integral over -1 to 1 is in fact done from 0 to 1 with a cosh, then doubled to even things out
            for (unsigned int i = 0;i < integrand.size();++i)
            {
                double t = i*step;

                integrand[i] = rNu*kappa*t*t/2.0
                        + anima::log_bessel_i(0, rNu*(kappa + 1.0)*wperp*std::sqrt(1-t*t) / sphereRadius)
                        - kappa
                        - ((kappa + 1.0)*wperpSq + wparaSq) / gaussDenom
                        - rNu*(kappa + 1.0)/2.0;

                double internalSinhVal = (kappa + rNu*wpara / sphereRadius) * t;
                if (std::abs(internalSinhVal) > 50)
                    integrand[i] += std::abs(internalSinhVal) - std::log(2.0);
                else
                    integrand[i] += std::log(std::cosh(internalSinhVal));

                integrand[i] = std::exp(integrand[i]);
            }

            double integral = 0;
            for (unsigned int i = 0;i < integrand.size() - 1;++i)
                integral += step * (integrand[i] + integrand[i+1])/2.0;

            resVal = 2.0 * C * integral;
        }
    }

    return resVal;
}

template <class VectorType, class ScalarType> double ComputeSymmetricPDF(const VectorType &w, const VectorType &mu,
                                                                         const ScalarType &kappa, const ScalarType &d, const ScalarType &nu,
                                                                         const ScalarType &step, std::vector<ScalarType> &integrand)
{
    double invNu = 1.0 - nu;

    if (invNu < 0.0)
        invNu = 0.0;

    if (invNu > 1.0)
        invNu = 1.0;

    if (invNu < 1e-4 || d < 1e-16)
        return 0;

    double rNu = nu / invNu;
    double gaussDenom = 2.0*invNu*d;
    double C = (kappa + 1.0) / pow(gaussDenom*M_PI,1.5);
    double wNorm = anima::ComputeNorm(w);
    double sphereRadius = std::sqrt(nu * d);

    double resVal;

    if (nu < 1e-8)
    {
        double wpara = anima::ComputeScalarProduct(mu, w);
        double wparaSq = wpara * wpara;
        double wperpSq = wNorm * wNorm - wparaSq;
        resVal = C * exp(-((kappa + 1.0)*wperpSq+wparaSq)/gaussDenom);
    }
    else
    {
        if (kappa < 1e-4)
        {
            wNorm /= sphereRadius;

            C *= ((1.0 - std::exp(-2.0*rNu*wNorm)) / (2.0*rNu*wNorm));
            resVal = C * std::exp(-rNu*wNorm*wNorm/2.0
                                  +rNu*wNorm
                                  -rNu/2.0);
        }
        else
        {
            C *= (kappa / (1.0 - exp(-2.0*kappa)));

            double wpara = anima::ComputeScalarProduct(mu, w);
            double wparaSq = wpara * wpara;
            double wperpSq = wNorm * wNorm - wparaSq;

            if (wperpSq < 0)
                wperpSq = 0;

            double wperp = std::sqrt(wperpSq);

            // Integral over -1 to 1 is in fact done from 0 to 1 with a cosh, then doubled to even things out
            // That plus symmetrization leads to two cosh values
            for (unsigned int i = 0;i < integrand.size();++i)
            {
                double t = i*step;

                integrand[i] = rNu*kappa*t*t/2.0
                        + anima::log_bessel_i(0, rNu*(kappa + 1.0)*wperp*std::sqrt(1-t*t) / sphereRadius)
                        - kappa
                        - ((kappa + 1.0)*wperpSq + wparaSq) / gaussDenom
                        - rNu*(kappa + 1.0)/2.0;

                if (std::abs(rNu*wpara*t / sphereRadius) > 50)
                    integrand[i] += std::abs(rNu*wpara*t / sphereRadius) - std::log(2.0);
                else
                    integrand[i] += std::log(std::cosh(rNu*wpara*t / sphereRadius));

                if (std::abs(kappa*t) > 50)
                    integrand[i] += std::abs(kappa*t) - std::log(2.0);
                else
                    integrand[i] += std::log(std::cosh(kappa*t));

                integrand[i] = std::exp(integrand[i]);
            }

            double integral = 0;
            for (unsigned int i = 0;i < integrand.size() - 1;++i)
                integral += step * (integrand[i] + integrand[i+1])/2.0;

            resVal = 2.0 * C * integral;
        }
    }

    return resVal;
}

template <class T1, class VectorType> double ComputeSymmetricCDF(const VectorType &direction, const T1 &kappa, const T1 &lambda, const T1 &nu,
                                                                 const T1 &bvalue, const VectorType &gradient)
{
    double resVal;

    double vmfLambda = nu * lambda;
    double gaussLambda = (1.0 - nu) * lambda;

    double muG = anima::ComputeScalarProduct(direction, gradient);

    resVal = std::exp( -bvalue * gaussLambda * (1 + kappa * muG * muG) / (kappa + 1) );

    double sqrt2bLambda = std::sqrt(2.0 * bvalue * vmfLambda);
    double kappaSq = kappa * kappa;

    double ReZ = kappaSq - sqrt2bLambda * sqrt2bLambda;
    double ImZ = 2.0 * kappa * sqrt2bLambda * muG;

    double condition = ReZ + std::sqrt(ReZ*ReZ + ImZ*ImZ);
    if (condition < 0)
        condition = 0;
    bool omega = condition < 1e-2;

    if (omega)
    {
        if (ReZ > 0)
            ReZ = 0;
        resVal *= (anima::SinOverId(sqrt(-ReZ)) / anima::ShOverId(kappa));
    }
    else
    {
        double alpha = std::sqrt(condition / 2.0);
        double beta = ImZ / std::sqrt(2.0 * condition);
        double tmp = anima::ShRatio(kappa, alpha, beta);
        resVal *= tmp;
    }

    return resVal;
}

template <class T1, class T2, class VectorType> double ComputeMixtureSymmetricCDF(const unsigned int &NbComponents,
                                                                                  const T2 *Directions, const T1 *Kappa,
                                                                                  const T1 *Lambda, const T1 *Nu, const T1 *W,
                                                                                  const T1 &b, const VectorType &grad)
{
    /* Isotropic component */
    double w0 = W[NbComponents];
    double lambda_iso = Lambda[NbComponents];
    double nu_iso = Nu[NbComponents];

    double resVal = 0;

    if (nu_iso == 0)
        resVal += w0 * std::exp(-b * lambda_iso);
    else
        resVal += w0 * std::exp(-b * (1.0 - nu_iso) * lambda_iso) * anima::SinOverId(sqrt(2.0 * b * nu_iso * lambda_iso));

    /* Other NbComponents component(s) of the mixture */
    for (unsigned int i = 0; i < NbComponents; ++i)
    {
        double w = W[i];

        if (w > 0)
        {
            double DDISingle = 0;

            if (Nu[i] == 0)
            {
                double cosAngle = anima::ComputeScalarProduct(grad,Directions[i]);
                double kVal = Kappa[i];
                DDISingle = std::exp( -b * Lambda[i] / (kVal + 1.0) * (1.0 + cosAngle * cosAngle * kVal) );
            }
            else
            {
                DDISingle = ComputeSymmetricCDF(Directions[i], Kappa[i], Lambda[i], Nu[i], b, grad);
            }

            resVal += w * DDISingle;
        }
    }

    return std::abs(resVal);
}

} // end namespace ddi_distribution
