#include <cmath>

#include "animaODFFunctions.h"

#include <iostream>
#include <fstream>

#include <animaVectorOperations.h>
#include <animaBaseTensorTools.h>
#include <itkSymmetricEigenAnalysis.h>

#include <boost/math/special_functions/factorials.hpp>

namespace anima
{

std::vector <std::vector <double> > InitializeSampleDirections(unsigned int nbTheta, unsigned int nbPhi,
                                                               std::string sampleDirFileName)
{
    std::vector <std::vector <double> > resVal;

    if (sampleDirFileName == "")
    {
        if ((nbTheta == 0)&&(nbPhi == 0))
            return resVal;

        std::vector <double> tmpVector(2,0);
        for (unsigned int i = 0;i < nbTheta;++i)
        {
            // Theta varies between 0 and pi/2
            tmpVector[0] = M_PI/(4*nbTheta) + i*M_PI/(2*nbTheta);
            for (unsigned int j = 0;j < nbPhi;++j)
            {
                // Phi varies between 0 and 2 pi
                tmpVector[1] = M_PI/nbPhi + 2*i*M_PI/nbPhi;
                resVal.push_back(tmpVector);
            }
        }
    }
    else
    {
        std::ifstream inputDirectionsFile(sampleDirFileName.c_str());
        std::string errorMsg = "Could not open sample directions file (";
        errorMsg += sampleDirFileName;
        errorMsg += ")";

        if (!inputDirectionsFile.is_open())
            throw itk::ExceptionObject(__FILE__, __LINE__,errorMsg,ITK_LOCATION);

        std::vector <double> gradTmp(3,0);
        std::vector <double> sphericalGrad(2,0);

        std::cout << "Using sample directions from file " << sampleDirFileName.c_str() << "..." << std::endl;

        while (!inputDirectionsFile.eof())
        {
            char tmpStr[2048];
            inputDirectionsFile.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            sscanf(tmpStr,"%lf %lf %lf",&gradTmp[0],&gradTmp[1],&gradTmp[2]);
            anima::TransformCartesianToSphericalCoordinates(gradTmp,sphericalGrad);
            resVal.push_back(sphericalGrad);
        }

        inputDirectionsFile.close();
    }

    return resVal;
}

void GetEulerAnglesFromRotationMatrix(vnl_matrix <double> &R, std::vector <double> &resVal)
{
    // This is not really clean. There are still undetermined cases.
    // R is actually R^T so adapting the values given by Geng et al.
    // TO DO : or is it ???
    resVal.resize(3);

    if (R(2,2) != 1)
    {
        resVal[0] = std::atan2(R(1,2),R(0,2)); // alpha
        resVal[1] = std::acos(R(2,2)); // beta
        resVal[2] = std::atan2(R(2,1),-R(2,0)); // gamma
    }
    else
    {
        if (R(1,0) < 0)
            resVal[0] = - std::acos(R(0,0));
        else
            resVal[0] = std::acos(R(0,0));

        resVal[1] = 0;
        resVal[2] = 0;
    }
}

double GetDValue(unsigned int l, int m, int mp, double angle)
{
    int globalSign = 1;
    if ((m - mp) % 2 == 1)
        globalSign = -1;

    long double factor = sqrt(boost::math::factorial<double>(l+m) * boost::math::factorial<double>(l-m) *
                              boost::math::factorial<double>(l+mp) * boost::math::factorial<double>(l-mp));

    int lowBoundK = 0;
    if (lowBoundK < mp - m)
        lowBoundK = mp - m;

    int supBoundK = l + mp;
    if (supBoundK > (int)l - m)
        supBoundK = (int)l - m;

    long double resVal = 0;

    for (int k = lowBoundK;k <= supBoundK;++k)
    {
        long double cosVal = 1, sinVal = 1;

        long double factDenom = boost::math::factorial<double>(k) * boost::math::factorial<double>(l - m - k) *
        boost::math::factorial<double>(l + mp - k) * boost::math::factorial<double>(m - mp + k);

        int powNum = 2*(l-k) - m + mp;

        double cosValTmp = cos(angle/2);
        for (int i = 1;i <= abs(powNum);++i)
            cosVal *= cosValTmp;

        if (powNum < 0)
            cosVal = 1.0/cosVal;

        powNum = 2*k + m - mp;

        long double sinValTmp = sin(angle/2);
        for (int i = 1;i <= abs(powNum);++i)
            sinVal *= sinValTmp;

        if (powNum < 0)
            sinVal = 1.0/sinVal;

        int kSign = 1;
        if (k % 2 == 1)
            kSign = -1;

        resVal += kSign*cosVal*sinVal/factDenom;
    }

    resVal *= factor*globalSign;
    return resVal;
}

void EstimateLocalODFRotationMatrix(vnl_matrix <double> &resVal, unsigned int l,
                                    double alpha, double beta, double gamma)
{
    unsigned int sizeMat = 2*l + 1;
    resVal.set_size(sizeMat,sizeMat);

    // Precompute d-values for use in inner loops
    vnl_matrix <double> dValues(2*l+1,2*l+1);

    dValues(l,l) = GetDValue(l,0,0,beta);
    for (unsigned int m = 1;m <= l;++m)
        dValues(l+m,l) = GetDValue(l,m,0,beta);

    for (unsigned int mp = 1;mp <= l;++mp)
    {
        for (unsigned int m = 1;m <= l;++m)
        {
            dValues(m+l,mp+l) = GetDValue(l,m,mp,beta);
            dValues(l-m,mp+l) = GetDValue(l,-m,mp,beta);
        }
    }

    for (int m = -l;m < 0;++m)
    {
        unsigned int absm = abs(m);

        int signm = 1;
        if (m % 2 == 1)
            signm = -1;

        for (int mp = -l;mp < 0;++mp)
        {
            int signmp = 1;
            if (mp % 2 == 1)
                signmp = -1;

            // Case all inferior to 0 ->double checked
            unsigned int absmp = abs(mp);
            resVal(m+l,mp+l) = dValues(absmp+l,absm+l)*cos(m*alpha+mp*gamma) + signmp*dValues(m+l,absmp+l)*cos(m*alpha-mp*gamma);
        }

        for (int mp = 1;mp <= (int)l;++mp)
        {
            // Case m < 0, mp > 0
            resVal(m+l,mp+l) = dValues(m+l,mp+l)*sin(m*alpha+mp*gamma) - signm*dValues(l+absm,l+mp)*sin(m*alpha-mp*gamma);
        }
    }

    for (int m = 1;m <= (int)l;++m)
    {
        int signm1 = 1;
        if (m % 2 == 0)
            signm1 = -1;

        for (int mp = -l;mp < 0;++mp)
        {
            int signall = 1;
            if ((m+mp) % 2 == 1)
                signall = -1;

            int signmp1 = 1;
            if (mp % 2 == 0)
                signmp1 = -1;

            // Case m > 0, mp < 0
            resVal(m+l,mp+l) = - signall*dValues(l-m,l-mp)*sin(m*alpha+mp*gamma) + signmp1*dValues(l+m,l-mp)*sin(m*alpha-mp*gamma);
        }


        for (int mp = 1;mp <= (int)l;++mp)
        {
            // Case m > 0, mp > 0 ->double checked
            resVal(m+l,mp+l) = dValues(l+m,l+mp)*cos(m*alpha+mp*gamma) + signm1*dValues(l-m,l+mp)*cos(m*alpha-mp*gamma);
        }
    }

    // This leaves only cases where one or all of m or mp are 0

    double sqrt2 = sqrt(2.0);
    for (int m = 1;m <= (int)l;++m)
    {
        double dValue = dValues(l+m,l);
        double oppDValue = dValue;
        if (m % 2 == 1)
            oppDValue = -dValue;

        resVal(l-m,l) = sqrt2*oppDValue*cos(m*alpha);
        resVal(l,l-m) = sqrt2*dValue*cos(m*gamma);

        resVal(m+l,l) = - sqrt2*dValue*sin(m*alpha);
        resVal(l,m+l) = sqrt2*oppDValue*sin(m*gamma);
    }

    resVal(l,l) = GetDValue(l,0,0,beta);
}

} // end of namespace anima
