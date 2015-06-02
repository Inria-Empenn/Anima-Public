#pragma once

#include "animaGradientFileReader.h"

#include <fstream>
#include <sstream>
#include <itkExceptionObject.h>

#include <animaVectorOperations.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace anima
{

    template <class GradientType, class BValueScalarType>
    GradientFileReader<GradientType,BValueScalarType>::
    GradientFileReader()
    {
        m_GradientFileName = "";
        m_BValueBaseString = "";

        m_B0ValueThreshold = 10;
        m_TotalNumberOfDirections = 0;

        m_GradientIndependentNormalization = true;
        m_Modified = false;
    }

    template <class GradientType, class BValueScalarType>
    std::vector <bool>
    GradientFileReader<GradientType,BValueScalarType>::
    GetB0Directions()
    {
        if (m_Modified)
            this->Update();

        std::vector <bool> resVal(m_TotalNumberOfDirections,false);

        for (unsigned int i = 0;i < m_TotalNumberOfDirections;++i)
        {
            if (m_BValues[i] <= m_B0ValueThreshold)
                resVal[i] = true;
        }

        return resVal;
    }

    template <class GradientType, class BValueScalarType>
    typename GradientFileReader<GradientType,BValueScalarType>::GradientVectorType &
    GradientFileReader<GradientType,BValueScalarType>::
    GetGradients()
    {
        if (m_Modified)
            this->Update();

        return m_Gradients;
    }

    template <class GradientType, class BValueScalarType>
    typename GradientFileReader<GradientType,BValueScalarType>::BValueVectorType &
    GradientFileReader<GradientType,BValueScalarType>::
    GetBValues()
    {
        if (m_Modified)
            this->Update();

        return m_BValues;
    }

    template <class GradientType, class BValueScalarType>
    unsigned int
    GradientFileReader<GradientType,BValueScalarType>::
    GetTotalNumberOfDirections()
    {
        if (m_Modified)
            this->Update();

        return m_TotalNumberOfDirections;
    }

    template <class GradientType, class BValueScalarType>
    void
    GradientFileReader<GradientType,BValueScalarType>::
    Update()
    {
        if (!m_Modified)
            return;

        if (boost::algorithm::ends_with(m_GradientFileName, "bvecs") || boost::algorithm::ends_with(m_GradientFileName, "bvec"))
            this->ReadGradientsFromBVecFile();
        else
            this->ReadGradientsFromTextFile();

        if (m_BValueBaseString != "")
        {
            this->ReadBValues();

            // Normalize gradient directions
            for (unsigned int i = 0;i < m_TotalNumberOfDirections;++i)
            {
                if ((m_BValues[i] <= m_B0ValueThreshold)||(anima::ComputeNorm(m_Gradients[i]) == 0))
                {
                    m_BValues[i] = 0;

                    for (unsigned int j = 0;j < 3;++j)
                        m_Gradients[i][j] = 0;

                    continue;
                }

                double norm = 1;
                if (!m_GradientIndependentNormalization)
                    norm = anima::ComputeNorm(m_Gradients[i]);

                anima::Normalize(m_Gradients[i],m_Gradients[i]);

                if (!m_GradientIndependentNormalization)
                    m_BValues[i] *= norm;
            }
        }

        m_Modified = false;
    }

    template <class GradientType, class BValueScalarType>
    void
    GradientFileReader<GradientType,BValueScalarType>::
    ReadGradientsFromBVecFile()
    {
        std::ifstream gradFile(m_GradientFileName.c_str());
        m_TotalNumberOfDirections = 0;

        if (!gradFile.is_open())
            throw itk::ExceptionObject(__FILE__, __LINE__,"Could not open gradient file as a bvec file",ITK_LOCATION);

        m_Gradients.clear();
        std::vector < std::vector <double> > trVecs(3);

        for (unsigned int i = 0;i < 3;++i)
        {
            if (gradFile.eof())
                throw itk::ExceptionObject(__FILE__, __LINE__,"BVec file is missing data",ITK_LOCATION);

            char tmpStr[8192];
            tmpStr[0] = '\0';

            while((strcmp(tmpStr,"") == 0)&&(!gradFile.eof()))
                gradFile.getline(tmpStr,8192);

            std::string workStr(tmpStr);
            boost::algorithm::trim_right(workStr);

            std::istringstream iss(workStr);
            std::string shortStr;
            do
            {
                iss >> shortStr;
                trVecs[i].push_back(boost::lexical_cast<BValueScalarType> (shortStr));
            }
            while (!iss.eof());
        }

        m_TotalNumberOfDirections = trVecs[0].size();

        m_Gradients.resize(m_TotalNumberOfDirections);
        for (unsigned int i = 0;i < m_TotalNumberOfDirections;++i)
        {
            this->InitializeEmptyGradient(m_Gradients[i]);
            for (unsigned int j = 0;j < 3;++j)
                m_Gradients[i][j] = trVecs[j][i];
        }
    }

    template <class GradientType, class BValueScalarType>
    void
    GradientFileReader<GradientType,BValueScalarType>::
    ReadGradientsFromTextFile()
    {
        std::ifstream gradFile(m_GradientFileName.c_str());
        m_TotalNumberOfDirections = 0;

        if (!gradFile.is_open())
        {
            throw itk::ExceptionObject(__FILE__, __LINE__,"Could not open gradient file as a text file",ITK_LOCATION);
        }

        GradientType gradTmp;
        this->InitializeEmptyGradient(gradTmp);
        m_Gradients.clear();

        while (!gradFile.eof())
        {
            char tmpStr[2048];
            gradFile.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            float a0, a1, a2;
            sscanf(tmpStr,"%f %f %f",&a0,&a1,&a2);
            gradTmp[0] = a0;
            gradTmp[1] = a1;
            gradTmp[2] = a2;

            m_Gradients.push_back(gradTmp);
        }

        m_TotalNumberOfDirections = m_Gradients.size();

        gradFile.close();
    }

    template <class GradientType, class BValueScalarType>
    void
    GradientFileReader<GradientType,BValueScalarType>::
    ReadBValues()
    {
        m_BValues.resize(m_TotalNumberOfDirections);

        std::ifstream bvalFile(m_BValueBaseString.c_str());

        if (!bvalFile.is_open())
        {
            for (unsigned int i = 0;i < m_TotalNumberOfDirections;++i)
                m_BValues[i] = boost::lexical_cast<BValueScalarType>(m_BValueBaseString);

            return;
        }

        if (boost::algorithm::ends_with(m_BValueBaseString, "bvals") || boost::algorithm::ends_with(m_BValueBaseString, "bval")) // bval file
        {
            char tmpStr[8192];
            bvalFile.getline(tmpStr,8192);
            std::string workStr(tmpStr);
            boost::algorithm::trim_right(workStr);

            std::istringstream iss(workStr);
            std::string shortStr;
            unsigned int i = 0;
            do
            {
                iss >> shortStr;
                m_BValues[i] = boost::lexical_cast<BValueScalarType> (shortStr);
                ++i;
            }
            while (!iss.eof());
        }
        else // text file
        {
            unsigned int i = 0;
            while (!bvalFile.eof())
            {
                char tmpStr[2048];
                bvalFile.getline(tmpStr,2048);

                if (strcmp(tmpStr,"") == 0)
                    continue;

                m_BValues[i] = boost::lexical_cast<BValueScalarType> (tmpStr);
                ++i;
            }
        }

        bvalFile.close();
    }

} // end namespace anima
