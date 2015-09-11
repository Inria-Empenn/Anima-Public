#include <animaRefactoredBaseBMRegistrationMethod.h>
#include <itkVectorImage.h>

int main(int argc, char **argv)
{
    typedef anima::RefactoredBaseBMRegistrationMethod <itk::Image <float, 3> > RegistrationType;

    RegistrationType::Pointer tmpPtr = RegistrationType::New();

    tmpPtr->Update();

    typedef anima::RefactoredBaseBMRegistrationMethod <itk::VectorImage <float, 3> > VectorRegistrationType;

    VectorRegistrationType::Pointer tmpPtr2 = VectorRegistrationType::New();

    tmpPtr2->Update();

    return 0;
}
