#pragma once

#include <vnl/vnl_vector_fixed.h>
#include <itkPoint.h>
#include <itkVector.h>
#include <vector>
#include <string>

namespace anima
{
    
    template <class GradientType = std::vector <double>, class BValueScalarType = double>
    class GradientFileReader
    {
    public:
//        typedef std::vector <GradientScalarType> GradientType;
        typedef typename std::vector <GradientType> GradientVectorType;
        typedef std::vector <BValueScalarType> BValueVectorType;
       
        GradientFileReader();
        ~GradientFileReader() {}
        
        void SetGradientFileName(std::string fName) {m_GradientFileName = fName;this->SetModified();}
        void SetBValueBaseString(std::string base) {m_BValueBaseString = base;this->SetModified();}
        void SetB0ValueThreshold(double thr) {m_B0ValueThreshold = thr;this->SetModified();}
        void SetGradientIndependentNormalization(bool normGrads) {m_GradientIndependentNormalization = normGrads;this->SetModified();}
        void SetModified() {m_Modified = true;}
        
        std::vector <bool> GetB0Directions();
        GradientVectorType &GetGradients();
        BValueVectorType &GetBValues();
        unsigned int GetTotalNumberOfDirections();
        
        void Update();
        
    protected:
        template <class T> void InitializeEmptyGradient(std::vector <T> &vec)
        {
            vec.resize(3);
            std::fill(vec.begin(),vec.end(),0);
        }
        
        template <class T> void InitializeEmptyGradient(itk::Point <T,3> &vec)
        {
            vec.Fill(0);
        }
        
        template <class T> void InitializeEmptyGradient(vnl_vector_fixed <T,3> &vec)
        {
            vec.fill(0);
        }
        
        template <class T> void InitializeEmptyGradient(itk::Vector <T,3> &vec)
        {
            vec.Fill(0);
        }
        
    private:
        void ReadGradientsFromBVecFile();
        void ReadGradientsFromTextFile();
        void ReadBValues();
        
        GradientVectorType m_Gradients;
        BValueVectorType m_BValues;
        
        std::string m_GradientFileName;
        std::string m_BValueBaseString;
        
        unsigned int m_TotalNumberOfDirections;
        
        double m_B0ValueThreshold;
        // If set to true, gradients are normalized to 1 without changing b-values, otherwise yes
        bool m_GradientIndependentNormalization;
        
        // Internal value to trigger update
        bool m_Modified;
    };
    
} // end namespace anima

#include "animaGradientFileReader.hxx"
