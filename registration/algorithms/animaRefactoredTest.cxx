#include <animaRefactoredAsymmetricBMRegistrationMethod.h>
#include <itkVectorImage.h>

int main(int argc, char **argv)
{
    typedef anima::RefactoredAsymmetricBMRegistrationMethod <itk::Image <float, 3> > RegistrationType;

    RegistrationType::Pointer tmpPtr = RegistrationType::New();

    tmpPtr->Update();

    typedef anima::RefactoredAsymmetricBMRegistrationMethod <itk::VectorImage <float, 3> > VectorRegistrationType;

    VectorRegistrationType::Pointer tmpPtr2 = VectorRegistrationType::New();

    tmpPtr2->Update();

    return 0;
}
