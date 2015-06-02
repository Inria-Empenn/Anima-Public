#pragma once

#include <itkGeneralTransform.h>

namespace anima
{

    template <class TScalarType=float, unsigned int NDimensions=3>
    class TransformSeriesReader
    {
    public:
        enum TransformationType
        {
            LINEAR,
            SVF_FIELD,
            DENSE_FIELD
        };

        struct TransformInformation
        {
            std::string fileName;
            TransformationType trType;
            bool invert;
        };

        typedef itk::GeneralTransform <TScalarType,NDimensions> OutputTransformType;
        typedef typename OutputTransformType::Pointer OutputTransformPointer;

        TransformSeriesReader();
        ~TransformSeriesReader();

        void SetInput(std::string const& trListName) {m_Input = trListName;}
        void SetInvertTransform(bool val) {m_InvertTransform = val;}

        void Update();

        OutputTransformType *GetOutputTransform() {return m_OutputTransform;}

    protected:
        void addLinearTransformation(std::string &fileName, bool invert);
        void addSVFTransformation(std::string &fileName, bool invert);
        void addDenseTransformation(std::string &fileName, bool invert);

    private:
        OutputTransformPointer m_OutputTransform;
        bool m_InvertTransform;

        std::string m_Input;
    };

} // end namespace itk

#include "animaTransformSeriesReader.hxx"
