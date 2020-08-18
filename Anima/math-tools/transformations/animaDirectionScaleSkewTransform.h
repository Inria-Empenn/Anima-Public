#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima
{
    
    template < class TScalarType=double >    // Data type for scalars (double or double)
    class DirectionScaleSkewTransform :
    public itk::MatrixOffsetTransformBase< TScalarType, 3 >
    {
    public:
        /** Standard class typedefs. */
        typedef DirectionScaleSkewTransform Self;
        typedef itk::MatrixOffsetTransformBase< TScalarType, 3 > Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;
        
        /** New macro for creation of through a Smart Pointer. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(DirectionScaleSkewTransform, MatrixOffsetTransformBase);
        
        /** Dimension of the space. */
        itkStaticConstMacro(SpaceDimension, unsigned int, 3);
        itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
        itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
        itkStaticConstMacro(ParametersDimension, unsigned int, 4);
        
        /** Parameters Type   */
        typedef typename Superclass::ParametersType         ParametersType;
        typedef typename Superclass::JacobianType           JacobianType;
        typedef typename Superclass::ScalarType             ScalarType;
        typedef typename Superclass::InputPointType         InputPointType;
        typedef typename Superclass::OutputPointType        OutputPointType;
        typedef typename Superclass::InputVectorType        InputVectorType;
        typedef typename Superclass::OutputVectorType       OutputVectorType;
        typedef typename Superclass::InputVnlVectorType     InputVnlVectorType;
        typedef typename Superclass::OutputVnlVectorType    OutputVnlVectorType;
        typedef typename Superclass::InputCovariantVectorType
        InputCovariantVectorType;
        typedef typename Superclass::OutputCovariantVectorType
        OutputCovariantVectorType;
        typedef typename Superclass::MatrixType             MatrixType;
        typedef itk::Matrix <ScalarType,4,4>                HomogeneousMatrixType;
        typedef typename Superclass::InverseMatrixType      InverseMatrixType;
        typedef typename Superclass::CenterType             CenterType;
        typedef typename Superclass::OffsetType             OffsetType;
        typedef typename Superclass::TranslationType        TranslationType;
        
        virtual void SetParameters(const ParametersType & parameters) override;
        virtual const ParametersType& GetParameters() const override;
        
        virtual void SetIdentity() override;
        
        virtual void SetGeometry(HomogeneousMatrixType &matrix, bool update = true);
        
        void SetLogScale(ScalarType scale, bool update = true);
        itkGetConstReferenceMacro(LogScale, ScalarType);

        void SetLogFirstSkew(ScalarType skew, bool update = true);
        itkGetConstReferenceMacro(LogFirstSkew, ScalarType);

        void SetLogSecondSkew(ScalarType skew, bool update = true);
        itkGetConstReferenceMacro(LogSecondSkew, ScalarType);

        void SetUniqueTranslation(ScalarType translation, bool update = true);
        itkGetConstReferenceMacro(UniqueTranslation, ScalarType);

        void SetUniqueDirection(unsigned int direction);
        itkGetConstReferenceMacro(UniqueDirection, unsigned int);
        
        virtual const itk::Vector <TScalarType,12>& GetLogVector() const {return m_LogVector;}
        virtual const vnl_matrix <TScalarType>& GetLogTransform() const {return m_LogTransform;}

    protected:
        DirectionScaleSkewTransform()
        : Superclass(ParametersDimension)
        {
            this->InitializeParameters();
        }

        DirectionScaleSkewTransform (unsigned int ParamersDim)
        : Superclass(ParamersDim)
        {
            this->InitializeParameters();
        }
        
        ~DirectionScaleSkewTransform() {}
        
        void PrintSelf(std::ostream &os, itk::Indent indent) const override;
        
        virtual void ComputeMatrix() override;
        virtual void GetInternalLogarithm();
        virtual void GetInternalExponential();
        virtual void InternalApplyGeometry(MatrixType &linearMatrix, OffsetType &offset);

        vnl_matrix <TScalarType> m_LogTransform, m_ExpTransform;
        HomogeneousMatrixType m_Geometry, m_GeometryInv;

    private:
        DirectionScaleSkewTransform(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented

        void InitializeParameters()
        {
            m_LogVector.Fill(itk::NumericTraits<TScalarType>::Zero);
            m_LogTransform.set_size(4,4);
            m_ExpTransform.set_size(4,4);

            m_LogScale = 0;
            m_LogFirstSkew = 0;
            m_LogSecondSkew = 0;
            m_UniqueTranslation = 0;
            m_UniqueDirection = 1; // 0 is the X axis, 1 the Y axis, 2 the Z axis of the image

            m_Geometry.SetIdentity();
            m_GeometryInv.SetIdentity();

            this->ComputeMatrix();
        }

        ScalarType  m_LogScale;
        ScalarType  m_LogFirstSkew, m_LogSecondSkew;
        ScalarType  m_UniqueTranslation;
        
        unsigned int m_UniqueDirection;
                
        itk::Vector <TScalarType,12> m_LogVector;
    }; // class DirectionScaleSkewTransform
    
    template <class TScalarType=double>    // Data type for scalars (double or double)
    class DirectionScaleTransform :
            public DirectionScaleSkewTransform <TScalarType>
    {
    public:
        /** Standard class typedefs. */
        typedef DirectionScaleTransform Self;
        typedef anima::DirectionScaleSkewTransform <TScalarType> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        typedef typename Superclass::ParametersType ParametersType;
        typedef typename Superclass::MatrixType MatrixType;
        typedef typename Superclass::OffsetType OffsetType;

        /** New macro for creation of through a Smart Pointer. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(DirectionScaleTransform, DirectionScaleSkewTransform);

        /** Dimension of the space. */
        itkStaticConstMacro(ParametersDimension, unsigned int, 2);

        virtual void SetParameters(const ParametersType & parameters) override;
        virtual const ParametersType& GetParameters() const override;

    protected:
        DirectionScaleTransform()
        : Superclass(ParametersDimension) {}

        virtual void GetInternalLogarithm() override;
        virtual void GetInternalExponential() override;
        virtual void InternalApplyGeometry(MatrixType &linearMatrix, OffsetType &offset) override;

    }; // class DirectionScaleTransform


    template <class TScalarType=double> // Data type for scalars (double or double)
    class DirectionTransform :
            public DirectionScaleSkewTransform <TScalarType>
    {
    public:
        /** Standard class typedefs. */
        typedef DirectionTransform Self;
        typedef anima::DirectionScaleSkewTransform <TScalarType> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        typedef typename Superclass::ParametersType ParametersType;
        typedef typename Superclass::MatrixType MatrixType;
        typedef typename Superclass::OffsetType OffsetType;

        /** New macro for creation of through a Smart Pointer. */
        itkNewMacro(Self)

        /** Run-time type information (and related methods). */
        itkTypeMacro(DirectionTransform, DirectionScaleSkewTransform)

        /** Dimension of the space. */
        itkStaticConstMacro(ParametersDimension, unsigned int, 1);

        virtual void SetParameters(const ParametersType & parameters) override;
        virtual const ParametersType& GetParameters() const override;

    protected:
        DirectionTransform()
        : Superclass(ParametersDimension) {}

        virtual void GetInternalLogarithm() override;
        virtual void GetInternalExponential() override;
        virtual void InternalApplyGeometry(MatrixType &linearMatrix, OffsetType &offset) override;

    }; // class DirectionTransform

}  // namespace anima

#include "animaDirectionScaleSkewTransform.hxx"
