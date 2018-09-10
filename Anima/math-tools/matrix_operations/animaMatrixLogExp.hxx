#pragma once

#include "animaMatrixLogExp.h"
#include <iostream>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_inverse.h>

namespace anima
{
template <class T> vnl_matrix <T> GetSquareRoot(const vnl_matrix <T> & m, const double precision, vnl_matrix <T> & resultM)
{
    // declarations

    double energy;
    unsigned int D = m.rows();

    vnl_matrix <T> Mk1(D,D), Mk(D,D), Yk1(D,D), Yk(D,D), interm(D,D), invMk(D,D), Id(D,D);
    int niter=1;
    int niterMax=60;

    // initializations

    Mk=m; Yk=m; Id=m;
    Id.set_identity();

    interm=Yk*Yk-m;
    energy=interm.frobenius_norm();

    //T n = m.Norm();
    T n = D;

    // loop

    while( (niter <= niterMax) && (energy > precision) )
    {
        T gamma=(T)fabs(pow(vnl_determinant(Mk),-1.0/(2.0*n)));
        vnl_matrix_inverse <T> inverter(Mk);
        invMk= inverter.inverse();
        Mk1=(Id+ (Mk*(gamma*gamma)+invMk*(1./(gamma*gamma))  )*0.5  )*0.5;
        Yk1=Yk*(Id+invMk*(1./(gamma*gamma)) )*0.5*gamma;

        Yk=Yk1; Mk = Mk1;

        niter++;

        interm=Yk*Yk-m;
        energy=interm.frobenius_norm();
    }

    if (niter == niterMax)
    {
        std::cout <<"\nWarning, max number of iteration reached in sqrt computation. Final energy is:" << energy << ".\n";
    }

    // std::cout << "niter="<<niter<<".\n";

    resultM=Mk;

    return Yk;
}

template <class T> vnl_matrix <T> GetPadeLogarithm(const vnl_matrix <T> & m,const int numApprox)
{
    unsigned int D = m.rows();

    vnl_matrix <T> log=m;
    vnl_matrix <T> Id(D,D);
    Id.set_identity();
    vnl_matrix <T> diff(D,D),interm2(D,D),interm3(D,D),sqr(D,D),cube(D,D);

    diff=Id-m;
    T energy=diff.frobenius_norm();

    if (energy > 0.5)
    {
        std::cout <<"Warning, matrix is not close enough to Id to call Pade approximation. Frobenius Distance=" << energy <<". Returning original matrix.\n";
        return log;
    }

    switch (numApprox)
    {

        case 1:
            interm2=diff*(-1.);
            interm3=(Id-diff*0.5);
            break;

        case 2:
            sqr=diff*diff;

            interm2=diff*(-1.0)+sqr*(0.5);
            interm3=(Id-diff+sqr);
            break;

        case 3:
            sqr=diff*diff;
            cube=sqr*diff;

            interm2=diff*(-1.0)+sqr+cube*(11.0/60.0);
            interm3=(Id-diff*(1.5)+sqr*0.6+cube*(-0.05));
            break;


        default:
            sqr=diff*diff;
            cube=sqr*diff;
            std::cout <<"\n\ndefault is taken in log Pade\n\n";
            interm2=diff*(-1.0)+sqr+cube*(11.0/60.0);
            interm3=(Id-diff*(1.5)+sqr*0.6+cube*(-0.05));
            break;
    }

    log=interm2 * vnl_inverse(interm3);

    return log;
}

/*
template <class T> vnl_matrix <T> GetLogarithm(const vnl_matrix <T> & m, const double precision, const int numApprox)
{
    // New version : this seems to work for non pathological matrices
    unsigned int ndim = m.rows();

    vnl_real_eigensystem eig(m);
    vnl_matrix < vcl_complex<T> > Vinv = vnl_matrix_inverse < vcl_complex<T> > ( eig.V );

    vnl_matrix < vcl_complex<T> > logMat = eig.D.asMatrix();
    for (unsigned int i = 0;i < ndim;++i)
        logMat(i,i) = vcl_log(logMat(i,i));

    logMat = eig.V * logMat * Vinv;

    vnl_matrix <T> resMat (ndim,ndim);
    for (unsigned int i = 0;i < ndim;++i)
        for (unsigned int j = 0;j < ndim;++j)
            resMat(i,j) = vcl_real <T> (logMat(i,j));

    return resMat;
}
*/

template <class T> vnl_matrix <T> GetLogarithm(const vnl_matrix <T> & m, const double precision, const int numApprox)
{
    unsigned int D = m.rows();
    //typedef double T;

    T factor=1.0;
    vnl_matrix <double> log=m;
    T energy;
    vnl_matrix <double> Id(D,D);
    Id.set_identity();
    vnl_matrix <double> interm(D,D),Yi(D,D),Yi1(D,D),resultM(D,D);
    vnl_matrix <double> matrix_sum=m;

    int niterMax=60; int niter=1;

    Yi=m;
    interm=Yi-Id;
    energy=interm.frobenius_norm();
    matrix_sum.fill(0);

    while ( (energy > 0.5) && (niter <= niterMax))
    {
        // std::cout <<"\nniter=" << niter <<", energy="<<energy<<".\n";

        Yi1=GetSquareRoot(Yi,precision,resultM);

        if ((!std::isfinite(resultM(0,0))) || (!std::isfinite(Yi(0,0))))
        {
            log.fill(0);
            return log;
        }

        matrix_sum+= (Id-resultM)*factor;

        Yi=Yi1;
        interm=Yi-Id;
        energy=interm.frobenius_norm();

        factor*=2.0;
        niter++;
    }

    log=GetPadeLogarithm(Yi,numApprox);
    if (!std::isfinite(log(0,0)))
    {
        log.fill(0);
        return log;
    }

    return log*factor+matrix_sum;
}


template <class T> vnl_matrix <T> GetExponential(const vnl_matrix <T> & m, const int numApprox)
{
    unsigned int D = m.rows();

    T factor;
    vnl_matrix <T> exp=m;
    vnl_matrix <T> Id(D,D); Id.set_identity();
    vnl_matrix <T> interm(D,D),interm2(D,D),interm3(D,D),sqr(D,D),cube(D,D);

    T norm=m.frobenius_norm();
    int k;

    if(norm > 1)
    {
        k=1 + (int) ceil(log(m.frobenius_norm())/log(2.0));
    }
    else if(norm >0.5)
    {
        k=1;
    }
    else
    {
        k=0;
    }

    //std::cout << "\n The famous k=" << ceil(log(m.GetVnlMatrix().frobenius_norm())/log(2.0)) << ". In effect norm=" << m.GetVnlMatrix().frobenius_norm()<< ".\n";

    factor=pow(2.0,(double) k);
    interm=m/factor;


    switch(numApprox)
    {
        case 1:
            interm2=Id+interm*0.5;
            interm3=Id-interm*0.5;
            break;

        case 2:
            sqr=interm*interm;

            interm2=Id+interm*0.5+sqr*(1.0/12.0);
            interm3=Id-interm*0.5+sqr*(1.0/12.0);
            break;

        case 3:
            sqr=interm*interm;
            cube=sqr*interm;

            interm2=Id+interm*0.5+sqr*(1.0/10.0)+cube*(1.0/120.0);
            interm3=Id-interm*0.5+sqr*(1.0/10.0)-cube*(1.0/120.0);;
            break;

        default:
            std::cout <<"\n\ndefault is taken in exp Pade\n\n";

            sqr=interm*interm;
            cube=sqr*interm;

            interm2=Id+interm*0.5+sqr*(1.0/10.0)+cube*(1.0/120.0);
            interm3=Id-interm*0.5+sqr*(1.0/10.0)-cube*(1.0/120.0);;
            break;
    }

    interm=interm2 * vnl_inverse(interm3);

    for(int i=1;i<=k;i++)
        interm*=interm;

    exp=interm;

    return exp;
}

template <class T> void Get3DRotationLogarithm(const vnl_matrix <T> &rotationMatrix, std::vector <T> &outputAngles)
{
    const unsigned int dataDim = 3;
    if (outputAngles.size() < dataDim)
        outputAngles.resize(dataDim);

    double rotationTrace = 0;
    for (unsigned int i = 0;i < dataDim;++i)
        rotationTrace += rotationMatrix(i,i);

    double cTheta = (rotationTrace - 1.0) / 2.0;
    double theta = std::acos(cTheta);
    double sTheta = std::sin(theta);

    unsigned int pos = 0;
    for (int i = dataDim - 1;i >= 0;--i)
    {
        for (int j = dataDim - 1;j > i;--j)
        {
            outputAngles[pos] = theta * (rotationMatrix(i,j) - rotationMatrix(j,i)) / (2.0 * sTheta);
            if (pos % 2 == 0)
                outputAngles[pos] *= -1;

            ++pos;
        }
    }
}

template <class T> void Get3DRotationExponential(const std::vector <T> &angles, vnl_matrix <T> &outputRotation)
{
    const unsigned int dataDim = 3;
    outputRotation.set_size(dataDim,dataDim);
    outputRotation.set_identity();

    double thetaValue = 0;
    for (unsigned int i = 0;i < dataDim;++i)
        thetaValue += angles[i] * angles[i];

    if (thetaValue == 0)
        return;

    double sqrtTheta = std::sqrt(thetaValue);
    double firstAddConstant = std::sin(sqrtTheta) / sqrtTheta;
    double secondAddConstant = (1.0 - std::cos(sqrtTheta)) / thetaValue;

    for (unsigned int i = 0;i < dataDim;++i)
    {
        for (unsigned int j = 0;j < dataDim;++j)
        {
            if (j != i)
            {
                outputRotation(i,i) -= secondAddConstant * angles[j] * angles[j];
                outputRotation(i,j) += secondAddConstant * angles[i] * angles[j];
            }
        }
    }

    // Add skew symmetric part
    unsigned int pos = dataDim - 1;
    for (unsigned int i = 0;i < dataDim;++i)
    {
        for (unsigned int j = i + 1;j < dataDim;++j)
        {
            if ((i + j) % 2 != 0)
            {
                outputRotation(i,j) -= angles[pos] * firstAddConstant;
                outputRotation(j,i) += angles[pos] * firstAddConstant;
            }
            else
            {
                outputRotation(i,j) += angles[pos] * firstAddConstant;
                outputRotation(j,i) -= angles[pos] * firstAddConstant;
            }

            --pos;
        }
    }
}

// Multi-threaded matrix logger
template <class TInputScalarType, class TOutputScalarType, unsigned int NDimensions, unsigned int NDegreesOfFreedom>
void
MatrixLoggerFilter<TInputScalarType,TOutputScalarType,NDimensions,NDegreesOfFreedom>::
Update()
{
    if (!m_Modified)
        return;

    m_Output.resize(m_InputTransforms.size());

    itk::MultiThreader::Pointer threaderLog = itk::MultiThreader::New();

    ThreadedLogData *tmpStr = new ThreadedLogData;
    tmpStr->loggerFilter = this;

    unsigned int actualNumberOfThreads = std::min(m_NumberOfThreads,(unsigned int)m_InputTransforms.size());

    threaderLog->SetNumberOfWorkUnits(actualNumberOfThreads);
    threaderLog->SetSingleMethod(this->ThreadedLogging,tmpStr);
    threaderLog->SingleMethodExecute();

    delete tmpStr;

    m_Modified = false;
}

template <class TInputScalarType, class TOutputScalarType, unsigned int NDimensions, unsigned int NDegreesOfFreedom>
ITK_THREAD_RETURN_TYPE
MatrixLoggerFilter<TInputScalarType,TOutputScalarType,NDimensions,NDegreesOfFreedom>::
ThreadedLogging(void *arg)
{
    itk::MultiThreader::WorkUnitInfo *threadArgs = (itk::MultiThreader::WorkUnitInfo *)arg;

    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int nbProcs = threadArgs->NumberOfWorkUnits;

    ThreadedLogData *tmpStr = (ThreadedLogData *)threadArgs->UserData;

    tmpStr->loggerFilter->InternalLogger(nbThread,nbProcs);

    return NULL;
}

template <class TInputScalarType, class TOutputScalarType, unsigned int NDimensions, unsigned int NDegreesOfFreedom>
void
MatrixLoggerFilter<TInputScalarType,TOutputScalarType,NDimensions,NDegreesOfFreedom>::
InternalLogger(unsigned int threadId, unsigned int numThreads)
{
    unsigned step = m_InputTransforms.size()/numThreads;

    unsigned startTrsf = threadId*step;
    unsigned endTrsf = (1+threadId)*step;

    if (threadId+1==numThreads) // Ensure we got up to the last trsf
        endTrsf = m_InputTransforms.size();

    typename BaseInputMatrixTransformType::MatrixType affinePart;
    itk::Vector <TInputScalarType,NDimensions> offsetPart;
    vnl_matrix <TInputScalarType> tmpMatrix(NDimensions+1,NDimensions+1,0), logMatrix(NDimensions+1,NDimensions+1,0);
    tmpMatrix(NDimensions,NDimensions) = 1;

    for (unsigned int i = startTrsf; i < endTrsf; ++i)
    {
        BaseInputMatrixTransformType *tmpTrsf = (BaseInputMatrixTransformType *)m_InputTransforms[i].GetPointer();
        affinePart = tmpTrsf->GetMatrix();
        offsetPart = tmpTrsf->GetOffset();

        for (unsigned int j = 0;j < NDimensions;++j)
        {
            tmpMatrix(j,NDimensions) = offsetPart[j];
            for (unsigned int k = 0;k < NDimensions;++k)
                tmpMatrix(j,k) = affinePart(j,k);
        }

        logMatrix = GetLogarithm(tmpMatrix);

        unsigned int pos = 0;
        if (m_UseRigidTransforms)
        {
            for (unsigned int j = 0;j < NDimensions;++j)
            {
                for (unsigned int k = j+1;k < NDimensions;++k)
                {
                    m_Output[i][pos] = logMatrix(j,k);
                    ++pos;
                }
            }

            for (unsigned int j = 0;j < NDimensions;++j)
            {
                m_Output[i][pos] = logMatrix(j,NDimensions);
                ++pos;
            }

        }
        else
        {
            for (unsigned int j = 0;j < NDimensions;++j)
            {
                for (unsigned int k = 0;k <= NDimensions;++k)
                {
                    m_Output[i][pos] = logMatrix(j,k);
                    ++pos;
                }
            }
        }
    }
}

} // end of namespace anima
