#pragma once

#include <itkVector.h>
#include <itkMatrix.h>
#include <itkSingleValuedNonLinearOptimizer.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/** \class BobyqaOptimizer
 * \brief BOBYQA Optimizer.
 *
 * \ingroup Numerics Optimizers
 */
class ANIMAOPTIMIZERS_EXPORT BobyqaOptimizer : public itk::SingleValuedNonLinearOptimizer
{
public:
    /** Standard class typedefs. */
    typedef BobyqaOptimizer           Self;
    typedef itk::Object                    Superclass;
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef itk::SingleValuedNonLinearOptimizer::ParametersType ParametersType;
    typedef itk::SingleValuedNonLinearOptimizer::ScalesType ScalesType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(BobyqaOptimizer, SingleValuedNonLinearOptimizer);


    /** Type of the Cost Function   */
    typedef  itk::SingleValuedCostFunction         CostFunctionType;
    typedef  CostFunctionType::Pointer        CostFunctionPointer;

    /** Set if the Optimizer should Maximize the metric */
    itkSetMacro( Maximize, bool );
    itkBooleanMacro( Maximize);
    itkGetConstReferenceMacro( Maximize, bool );

    /** Set/Get maximum iteration limit. */
    itkSetMacro( MaximumIteration, unsigned int );
    itkGetConstReferenceMacro( MaximumIteration, unsigned int );


    itkSetMacro(SpaceDimension, unsigned int);
    itkGetConstReferenceMacro(SpaceDimension, unsigned int);

    /** Set/Get Number of Sampling points   (between Parameters+2 and (Parameters+1)(Parameters+2)/2)*/
    itkSetMacro( NumberSamplingPoints, unsigned int );
    itkGetConstReferenceMacro( NumberSamplingPoints, unsigned int );


    /** Set/Get RhoBegin set the initial sampling of
     * parameter space */
    itkSetMacro( RhoBegin, double );
    itkGetConstReferenceMacro( RhoBegin, double );

    /** Set/Get RhoEnd set the final sampling distance of
     * parameter space */
    itkSetMacro( RhoEnd, double );
    itkGetConstReferenceMacro( RhoEnd, double );

    /** Return Current Value */
    itkGetConstReferenceMacro( CurrentCost, MeasureType );
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
    void StartOptimization();

    /** When users call StartOptimization, this value will be set false.
     * By calling StopOptimization, this flag will be set true, and
     * optimization will stop at the next iteration. */
    void StopOptimization()
    { m_Stop = true; }

    itkGetConstReferenceMacro(CatchGetValueException, bool);
    itkSetMacro(CatchGetValueException, bool);

    itkGetConstReferenceMacro(MetricWorstPossibleValue, double);
    itkSetMacro(MetricWorstPossibleValue, double);

    itkSetMacro(LowerBounds, ScalesType);
    itkSetMacro(UpperBounds, ScalesType);

    const std::string GetStopConditionDescription() const;

protected:
    BobyqaOptimizer();
    BobyqaOptimizer(const BobyqaOptimizer&);
    virtual ~BobyqaOptimizer();
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

    /** Set if the Metric should be maximized: Default = False */
    bool               m_Maximize;

    unsigned int       m_SpaceDimension;

    unsigned int m_MaximumIteration;
    unsigned int m_CurrentIteration;

    bool               m_CatchGetValueException;
    double             m_MetricWorstPossibleValue;

    unsigned int m_NumberSamplingPoints;
    double m_RhoBegin, m_RhoEnd;

    //        int  Parameters;
    BobyqaOptimizer::ParametersType px;
    ScalesType m_LowerBounds, m_UpperBounds;

    itkSetMacro(CurrentCost, double);

    /** Internal storage for the value type / used as a cache  */
    MeasureType        m_CurrentCost;

    /** this is user-settable flag to stop optimization.
     * when users call StartOptimization, this value will be set false.
     * By calling StopOptimization, this flag will be set true, and
     * optimization will stop at the next iteration. */
    bool               m_Stop;

    std::ostringstream m_StopConditionDescription;

protected: // From f2c...
    int optimize(long int &npt, double *x, double *xl, double *xu, double &rhobeg, double &rhoend,
                 long int &maxfun, double *w);
    int bobyqb_(long int &npt, double *x, double *xl, double *xu, double &rhobeg, double &rhoend,
                long int &maxfun, double *xbase, double *xpt, double *fval, double *xopt,
                double *gopt, double *hq, double *pq, double *bmat, double *zmat,
                long int &ndim, double *sl, double *su, double *xnew, double *xalt, double *d__,
                double *vlag, double *w);
    int altmov_(long int &npt, double *xpt, double *xopt, double *bmat, double *zmat, long int &ndim,
                double *sl, double *su, long int &kopt, long int &knew, double *adelt, double *xnew,
                double *xalt, double *alpha, double *cauchy, double *glag, double *hcol, double *w);
    int prelim_(long int &npt, double *x, double *xl, double *xu, double &rhobeg, long int &maxfun,
                double *xbase, double *xpt, double *fval, double *gopt, double *hq, double *pq,
                double *bmat, double *zmat, long int &ndim, double *sl, double *su, long int &nf,
                long int &kopt);
    int rescue_(long int &npt, double *xl, double *xu, long int &maxfun, double *xbase, double *xpt,
                double *fval, double *xopt, double *gopt, double *hq, double *pq, double *bmat,
                double *zmat, long int &ndim, double *sl, double *su, long int &nf, double *delta,
                long int &kopt, double *vlag, double *ptsaux, double *ptsid, double *w);
    int trsbox_(long int &npt, double *xpt, double *xopt, double *gopt, double *hq, double *pq,
                double *sl, double *su, double *delta, double *xnew, double *d__, double *gnew,
                double *xbdi, double *s, double *hs, double *hred, double *dsq, double *crvmin);
    int update_(long int &npt, double *bmat, double *zmat, long int &ndim, double *vlag,
                double *beta, double *denom, long int &knew, double *w);
};

} // end of namespace anima



