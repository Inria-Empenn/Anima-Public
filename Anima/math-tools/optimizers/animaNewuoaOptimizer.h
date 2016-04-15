/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkOptimizer.h,v $
 Language:  C++
 Date:      $Date: 2009-06-24 12:02:54 $
 Version:   $Revision: 1.39 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#pragma once

#include <itkVector.h>
#include <itkMatrix.h>
#include <itkSingleValuedNonLinearOptimizer.h>
#include "AnimaOptimizersExport.h"

namespace anima
{

/** \class NewuoaOptimizer
 * \brief Newuoa Optimizer.
 *
 * \ingroup Numerics Optimizers
 */
class ANIMAOPTIMIZERS_EXPORT NewuoaOptimizer : public itk::SingleValuedNonLinearOptimizer
{
public:
    /** Standard class typedefs. */
    typedef NewuoaOptimizer                Self;
    typedef itk::Object                    Superclass;
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef SingleValuedNonLinearOptimizer::ParametersType
    ParametersType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(NewuoaOptimizer, SingleValuedNonLinearOptimizer)


    /** Type of the Cost Function   */
    typedef  itk::SingleValuedCostFunction         CostFunctionType;
    typedef  CostFunctionType::Pointer        CostFunctionPointer;

    /** Set if the Optimizer should Maximize the metric */
    itkSetMacro(Maximize, bool)
    itkBooleanMacro(Maximize)
    itkGetConstReferenceMacro(Maximize, bool)

    /** Set/Get maximum iteration limit. */
    itkSetMacro(MaximumIteration, unsigned int)
    itkGetConstReferenceMacro(MaximumIteration, unsigned int)


    itkSetMacro(SpaceDimension, unsigned int)
    itkGetConstReferenceMacro(SpaceDimension, unsigned int)

    /** Set/Get Number of Sampling points   (between Parameters+2 and (Parameters+1)(Parameters+2)/2)*/
    itkSetMacro(NumberSamplingPoints, unsigned int)
    itkGetConstReferenceMacro(NumberSamplingPoints, unsigned int)

    /** Set/Get RhoBegin set the initial sampling of
     * parameter space */
    itkSetMacro(RhoBegin, double)
    itkGetConstReferenceMacro(RhoBegin, double)

    /** Set/Get RhoEnd set the final sampling distance of
     * parameter space */
    itkSetMacro(RhoEnd, double)
    itkGetConstReferenceMacro(RhoEnd, double)

    /** Return Current Value */
    itkGetConstReferenceMacro(CurrentCost, MeasureType)
    MeasureType GetValue() const
    {
        if (this->GetMaximize())
            return -this->GetCurrentCost();
        else
            return this->GetCurrentCost();
    }

    /** Return Current Iteration */
    itkGetConstReferenceMacro( CurrentIteration, unsigned int);

    /** Start optimization. */
    void StartOptimization() ITK_OVERRIDE;

    /** When users call StartOptimization, this value will be set false.
     * By calling StopOptimization, this flag will be set true, and
     * optimization will stop at the next iteration. */
    void StopOptimization() {m_Stop = true;}

    itkGetConstReferenceMacro(CatchGetValueException, bool)
    itkSetMacro(CatchGetValueException, bool)

    itkGetConstReferenceMacro(MetricWorstPossibleValue, double)
    itkSetMacro(MetricWorstPossibleValue, double)

    const std::string GetStopConditionDescription() const ITK_OVERRIDE;

protected:
    NewuoaOptimizer();
    NewuoaOptimizer(const NewuoaOptimizer&);
    virtual ~NewuoaOptimizer();
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    /** Set if the Metric should be maximized: Default = False */
    bool               m_Maximize;

    unsigned int       m_SpaceDimension;

    unsigned m_MaximumIteration;
    unsigned m_CurrentIteration;

    bool               m_CatchGetValueException;
    double             m_MetricWorstPossibleValue;

    unsigned m_NumberSamplingPoints;
    double m_RhoBegin, m_RhoEnd;

    //        int  Parameters;
    NewuoaOptimizer::ParametersType px;

    itkSetMacro(CurrentCost, double)

    /** Internal storage for the value type / used as a cache  */
    MeasureType        m_CurrentCost;

    /** this is user-settable flag to stop optimization.
     * when users call StartOptimization, this value will be set false.
     * By calling StopOptimization, this flag will be set true, and
     * optimization will stop at the next iteration. */
    bool               m_Stop;

    std::ostringstream m_StopConditionDescription;

protected: // From f2c...
    double optimize(double *w, long int npt, double *x,double rhobeg, double rhoend, long int maxfun, int &iter);
    double newuob_(long int npt, double *x, double rhobeg, double rhoend, long int maxfun, double *xbase, double *xopt, double *xnew,double *xpt, double *fval, double *gq, double *hq, double *pq, double *bmat, double *zmat, long int *ndim, double *d__, double *vlag, double *w, int &iter);
    int bigden_(long int npt, double *xopt,double *xpt, double *bmat, double *zmat, long int *idz,long int *ndim, long int *kopt, long int *knew, double *d__, double *w, double *vlag, double *beta, double *s, double *wvec, double *prod);
    int biglag_(long int npt, double *xopt, double *xpt, double *bmat, double *zmat, long int *idz, long int *ndim, long int *knew, double *delta, double *d__, double *alpha, double *hcol, double *gc, double *gd, double *s, double *w);
    int update_(long int npt, double *bmat,double *zmat, long int *idz, long int *ndim, double *vlag,double *beta, long int *knew, double *w);
    int trsapp_(long int npt, double *xopt,double *xpt, double *gq, double *hq, double *pq,double *delta, double *step, double *d__, double *g,double *hd, double *hs, double *crvmin);
};

} // end of namespace anima
