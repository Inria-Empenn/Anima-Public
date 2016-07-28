#pragma once

#include <itkVector.h>
#include <itkMatrix.h>
#include <itkSingleValuedNonLinearOptimizer.h>

#include <nlopt.h>
#include <AnimaOptimizersExport.h>

#include <animaNLOPTParametersConstraintFunction.h>

namespace anima
{
    /**
     * \section ITK Wrapper for the NLOPT library
     * Version 1.0.0\n
     *
     * This code was graciously provided to us by the CRL at Boston.
     *
     * Copyright (c) 2010-2011 Children's Hospital Boston.\n
     * Benoit Scherrer, CRL (Computational Radiology Laboratory), Harvard Medical School
     *
     * \section Example
     * The wrapper is very easy to use. Here is an example of how to minimimze a
     * function with MyCostFunction being a regular itk::SingleValuedCostFunction
     *
     * \code
     * 		// Create the cost function
     *		itk::MyCostFunction::Pointer costFun = itk::MyCostFunction::New();
     *
     *		// Create the optimizer
     *		itk::NLOPTOptimizers::Pointer optimizer = itk::NLOPTOptimizers::New();
     *		optimizer->SetCostFunction(costFun);
     *
     *		// Set the parameters
     *		optimizer->SetMaxEval(500);
     *		optimizer->SetAlgorithm(itk::NLOPTOptimizers::NLOPT_LN_BOBYQA);
     *
     *		// Set the initial position
     *		itk::NLOPTOptimizers::ParametersType pinit(2);
     *		pinit[0]=pinit[1] = 4;
     *		optimizer->SetInitialPosition(pinit);
     *
     *		optimizer->StartOptimization();
     *
     *		std::cout<<"Return code "<<(int)optimizer->GetErrorCode();
     *		std::cout<<" = "<<optimizer->GetErrorCodeDescription()<<endl;
     *		std::cout<<"Position = "<<optimizer->GetCurrentPosition()<<endl;
     *		std::cout<<"Value = "<<optimizer->GetCurrentCost()<<endl;
     *
     * \endcode
     *
     * For more information about NLOPT and the available algorithms: \n
     * http://ab-initio.mit.edu/wiki/index.php/NLopt
     *
     * \section License
     *
     * This software is licensed by the copyright holder under the terms of the
     * Open Software License version 3.0. \n
     * http://www.opensource.org/licenses/osl-3.0.php
     *
     * \section  Attribution Notice.
     *
     * This research was carried out in the Computational Radiology Laboratory of
     * Children's Hospital, Boston and Harvard Medical School.\n
     * http://www.crl.med.harvard.edu \n
     * For more information contact: simon.warfield@childrens.harvard.edu
     *
     * This research work was made possible by Grant Number R01 RR021885 (Principal
     * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
     * from the National Center for Research Resources (NCRR), a component of the
     * National Institutes of Health (NIH).
     *
     * \class	NLOPTOptimizers
     *
     * \brief	Implements an ITK wrapper for the NLOPT library.
     *
     * 			NLopt is a free/open-source library for nonlinear optimization, providing a common
     * 			interface for a number of different free optimization routines available online as
     * 			well as original implementations of various other algorithms.
     *
     * 			WEBSITE: http://ab-initio.mit.edu/wiki/index.php/NLopt \n
     * 			AUTHOR: Steven G. Johnson
     *
     *
     * 			This ITK wrapper was created by:\n Benoit Scherrer, CRL, Harvard Medical School. \n
     * 			Copyright (c) 2010-2011 Children's Hospital Boston
     *
     * 			\ingroup Numerics Optimizers
     *
     * \author	Benoit Scherrer
     * \date	July 2010
    *************************************************************************************************/
    class ANIMAOPTIMIZERS_EXPORT NLOPTOptimizers:
        public itk::SingleValuedNonLinearOptimizer
    {
    public:

        typedef NLOPTOptimizers                Self;			/**< Standard ITK "Self" typedef. */
        typedef SingleValuedNonLinearOptimizer Superclass;		/**< Standard ITK "Superclass" typedef. */
        typedef itk::SmartPointer<Self>             Pointer;			/**< Standard ITK "Pointer" typedef. */
        typedef itk::SmartPointer<const Self>       ConstPointer;	/**< Standard ITK "ConstPointer" typedef. */

        typedef SingleValuedNonLinearOptimizer::ParametersType ParametersType; /**< Parameter type */
        typedef anima::NLOPTParametersConstraintFunction ConstraintsFunctionType; /**< Base constraints (for inequality and equality) */

        /**********************************************************************************************//**
         * \enum	nlopt_algorithm
         *
         * \brief	Values that represent the nlopt_algorithm to use. Naming conventions:
         *
         * 			-  NLOPT_{G/L}{D/N}_ : global/local derivative/no-derivative optimization,
         * 			respectively
         *
         * 			-  _RAND : algorithms involve some randomization.
         *
         * 			-  *_NOSCAL : algorithms are *not* scaled to a unit hypercube
         * 			(i.e. they are sensitive to the units of x)
        *************************************************************************************************/
        typedef enum {
            NLOPT_GN_DIRECT = 0,
            NLOPT_GN_DIRECT_L,
            NLOPT_GN_DIRECT_L_RAND,
            NLOPT_GN_DIRECT_NOSCAL,
            NLOPT_GN_DIRECT_L_NOSCAL,
            NLOPT_GN_DIRECT_L_RAND_NOSCAL,
            
            NLOPT_GN_ORIG_DIRECT,
            NLOPT_GN_ORIG_DIRECT_L,
            
            NLOPT_GD_STOGO,
            NLOPT_GD_STOGO_RAND,
            
            NLOPT_LD_LBFGS_NOCEDAL,
            
            NLOPT_LD_LBFGS,
            
            NLOPT_LN_PRAXIS,
            
            NLOPT_LD_VAR1,
            NLOPT_LD_VAR2,
            
            NLOPT_LD_TNEWTON,
            NLOPT_LD_TNEWTON_RESTART,
            NLOPT_LD_TNEWTON_PRECOND,
            NLOPT_LD_TNEWTON_PRECOND_RESTART,
            
            NLOPT_GN_CRS2_LM,
            
            NLOPT_GN_MLSL,
            NLOPT_GD_MLSL,
            NLOPT_GN_MLSL_LDS,
            NLOPT_GD_MLSL_LDS,
            
            NLOPT_LD_MMA,
            
            NLOPT_LN_COBYLA,
            
            NLOPT_LN_NEWUOA,
            NLOPT_LN_NEWUOA_BOUND,
            
            NLOPT_LN_NELDERMEAD,
            NLOPT_LN_SBPLX,
            
            NLOPT_LN_AUGLAG,
            NLOPT_LD_AUGLAG,
            NLOPT_LN_AUGLAG_EQ,
            NLOPT_LD_AUGLAG_EQ,
            
            NLOPT_LN_BOBYQA,
            
            NLOPT_GN_ISRES,
            
            /* new variants that require local_optimizer to be set,
             not with older constants for backwards compatibility */
            NLOPT_AUGLAG,
            NLOPT_AUGLAG_EQ,
            NLOPT_G_MLSL,
            NLOPT_G_MLSL_LDS,
            
            NLOPT_LD_SLSQP,
            
            NLOPT_LD_CCSAQ,
            
            NLOPT_GN_ESCH
            
            //NLOPT_NUM_ALGORITHMS /* not an algorithm, just the number of them */
        } nlopt_algorithm;


        /**********************************************************************************************//**
         * \enum	nlopt_result
         *
         * \brief	Values that represent the result of a nlopt function.
        *************************************************************************************************/
        typedef enum {
            NLOPT_FAILURE = -1,				/**< generic failure code */
            NLOPT_INVALID_ARGS = -2,		/**< Invalid arguments (eg lower bounds bigger than upper bounds, unknown algorithm, etc) */
            NLOPT_OUT_OF_MEMORY = -3,		/**< Ran out of memory */
            NLOPT_ROUNDOFF_LIMITED = -4,	/**< Halted because roundoff errors limited progress */
            NLOPT_FORCED_STOP = -5,			/**< Halted because of a forced termination \sa StopOptimization */
            NLOPT_SUCCESS = 1,				/**< generic success code */
            NLOPT_STOPVAL_REACHED = 2,		/**< Optimization stopped because 'stopval' was reached */
            NLOPT_FTOL_REACHED = 3,			/**< Optimization stopped because 'ftol_rel' or 'ftol_abs' was reached */
            NLOPT_XTOL_REACHED = 4,			/**< Optimization stopped because 'xtol_rel' or 'xtol_abs' was reached */
            NLOPT_MAXEVAL_REACHED = 5,		/**< Optimization stopped because 'maxeval' was reached. */
            NLOPT_MAXTIME_REACHED = 6		/**< Optimization stopped because 'maxtime' was reached. */
        } nlopt_result;


        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(NLOPTOptimizers, SingleValuedNonLinearOptimizer );

        /** Type of the Cost Function   */
        typedef  itk::SingleValuedCostFunction         CostFunctionType;
        typedef  CostFunctionType::Pointer        CostFunctionPointer;

        /** NLOPT Algorithm to use */
        itkSetMacro( Algorithm, nlopt_algorithm );
        itkGetConstReferenceMacro( Algorithm, nlopt_algorithm );

        /** Returns the last error code of NLOPT */
        itkSetMacro( ErrorCode, nlopt_result );
        itkGetConstReferenceMacro( ErrorCode, nlopt_result );

        std::string GetErrorCodeDescription() const;

        bool	isSuccessful() const;


        /** Returns the current value */
        itkGetConstReferenceMacro( CurrentCost, MeasureType );
        /** Returns the current value */
        MeasureType GetValue() const { return this->GetCurrentCost(); }

        void SetLowerBoundParameters( const ParametersType& p );
        itkGetConstReferenceMacro( LowerBoundParameters, ParametersType );

        void SetUpperBoundParameters( const ParametersType& p );
        itkGetConstReferenceMacro( UpperBoundParameters, ParametersType );

        /** Set if the Optimizer should Maximize the metric */
        itkSetMacro( Maximize, bool );
        itkGetConstReferenceMacro( Maximize, bool );

        /** Stopping criteria:
          * Stop when an objective value of at least stopval is found: stop minimizing
          * when an objective value <= stopval is found, or stop maximizing a value >=
          * stopval is found. (Setting stopval to -HUGE_VAL for minimizing or +HUGE_VAL
          * for maximizing disables this stopping criterion.)
          *
          * \note Multiple stopping criteria for the optimization are
          * supported. The optimization halts whenever any one of these criteria is
          * satisfied. In some cases, the precise interpretation of the stopping
          * criterion depends on the optimization algorithm above (although NLopt
          * try to make them as consistent as reasonably possible), and some
          * algorithms do not support all of the stopping criteria.
          *
          * \note You do not need to use all of the stopping criteria! In most cases,
          * you only need one or two, and can omit the remainder (all criteria are
          * disabled by default).*/
        void SetStopVal(double stopval) { m_StopValSet=true; m_StopVal=stopval; }
        itkGetConstReferenceMacro( StopVal, double );

        /** Stopping criteria:
          * Set relative tolerance on function value: stop when an optimization step (or
          * an estimate of the optimum) changes the objective function value by less than
          * tol multiplied by the absolute value of the function value. (If there is any
          * chance that your optimum function value is close to zero, you might want to
          * set an absolute tolerance with nlopt_set_ftol_abs as well.) Criterion is
          * disabled if tol is non-positive.
          *
          * \note Multiple stopping criteria for the optimization are
          * supported. The optimization halts whenever any one of these criteria is
          * satisfied. In some cases, the precise interpretation of the stopping
          * criterion depends on the optimization algorithm above (although NLopt
          * try to make them as consistent as reasonably possible), and some
          * algorithms do not support all of the stopping criteria.
          *
          * \note You do not need to use all of the stopping criteria! In most cases,
          * you only need one or two, and can omit the remainder (all criteria are
          * disabled by default).*/
        itkSetMacro( FTolRel, double );
        itkGetConstReferenceMacro( FTolRel, double );

        /** Stopping criteria:
          * Set absolute tolerance on function value: stop when an optimization
          * step (or an estimate of the optimum) changes the function value by
          * less than tol. Criterion is disabled if tol is non-positive.
          *
          * \note Multiple stopping criteria for the optimization are
          * supported. The optimization halts whenever any one of these criteria is
          * satisfied. In some cases, the precise interpretation of the stopping
          * criterion depends on the optimization algorithm above (although NLopt
          * try to make them as consistent as reasonably possible), and some
          * algorithms do not support all of the stopping criteria.
          *
          * \note You do not need to use all of the stopping criteria! In most cases,
          * you only need one or two, and can omit the remainder (all criteria are
          * disabled by default).*/
        itkSetMacro( FTolAbs, double );
        itkGetConstReferenceMacro( FTolAbs, double );

        /** Stopping criteria:
          * Set relative tolerance on optimization parameters: stop when an optimization
          * step (or an estimate of the optimum) changes every parameter by less than
          * tol multiplied by the absolute value of the parameter. (If there is any
          * chance that an optimal parameter is close to zero, you might want to set
          * an absolute tolerance with nlopt_set_xtol_abs as well.) Criterion is
          * disabled if tol is non-positive.
          *
          * \note Multiple stopping criteria for the optimization are
          * supported. The optimization halts whenever any one of these criteria is
          * satisfied. In some cases, the precise interpretation of the stopping
          * criterion depends on the optimization algorithm above (although NLopt
          * try to make them as consistent as reasonably possible), and some
          * algorithms do not support all of the stopping criteria.
          *
          * \note You do not need to use all of the stopping criteria! In most cases,
          * you only need one or two, and can omit the remainder (all criteria are
          * disabled by default).*/
        itkSetMacro( XTolRel, double );
        itkGetConstReferenceMacro( XTolRel, double );

        /** Stopping criteria:
          * Set the absolute tolerance on optimization parameters.
          *
          * \todo Set the absolute tolerances for each parameter separately
          * by using nlopt_set_xtol_abs(nlopt_opt opt, const double* tol)
          *
          * \note Multiple stopping criteria for the optimization are
          * supported. The optimization halts whenever any one of these criteria is
          * satisfied. In some cases, the precise interpretation of the stopping
          * criterion depends on the optimization algorithm above (although NLopt
          * try to make them as consistent as reasonably possible), and some
          * algorithms do not support all of the stopping criteria.
          *
          * \note You do not need to use all of the stopping criteria! In most cases,
          * you only need one or two, and can omit the remainder (all criteria are
          * disabled by default).*/
        itkSetMacro( XTolAbs, double );
        itkGetConstReferenceMacro( XTolAbs, double );

        /** Stopping criteria:
          * Stop when the number of function evaluations exceeds maxeval. (This is not
          * a strict maximum: the number of function evaluations may exceed maxeval
          * slightly, depending upon the algorithm.) Criterion is disabled if maxeval
          * is non-positive.
          *
          * \note Multiple stopping criteria for the optimization are
          * supported. The optimization halts whenever any one of these criteria is
          * satisfied. In some cases, the precise interpretation of the stopping
          * criterion depends on the optimization algorithm above (although NLopt
          * try to make them as consistent as reasonably possible), and some
          * algorithms do not support all of the stopping criteria.
          *
          * \note You do not need to use all of the stopping criteria! In most cases,
          * you only need one or two, and can omit the remainder (all criteria are
          * disabled by default).*/
        itkSetMacro( MaxEval, int );
        itkGetConstReferenceMacro( MaxEval, int );

        /** Stopping criteria:
          * Stop when the optimization time (in seconds) exceeds maxtime. (This is not
          * a strict maximum: the time may exceed maxtime slightly, depending upon the
          * algorithm and on how slow your function evaluation is.) Criterion is disabled
          * if maxtime is non-positive
          *
          * \note Multiple stopping criteria for the optimization are
          * supported. The optimization halts whenever any one of these criteria is
          * satisfied. In some cases, the precise interpretation of the stopping
          * criterion depends on the optimization algorithm above (although NLopt
          * try to make them as consistent as reasonably possible), and some
          * algorithms do not support all of the stopping criteria.
          *
          * \note You do not need to use all of the stopping criteria! In most cases,
          * you only need one or two, and can omit the remainder (all criteria are
          * disabled by default).*/
        itkSetMacro( MaxTime, double );
        itkGetConstReferenceMacro( MaxTime, double );
        
        /** Vector storage size:
         * number of gradients to "remember" from previous optimization steps: increasing it increases the memory requirements but speeds convergence.*/
        itkSetMacro(VectorStorageSize, int);
        itkGetConstReferenceMacro(VectorStorageSize, int);
        
        /** Stochastic population size:
         * This is for random global optimization. It sets the number of random particles to propagate from.*/
        itkSetMacro(PopulationSize, int);
        itkGetConstReferenceMacro(PopulationSize, int);
        
        /** Set local optimizer for AUGLAG and MLSL algos */
        itkSetMacro( LocalOptimizer, nlopt_algorithm );
        itkGetConstReferenceMacro( LocalOptimizer, nlopt_algorithm );

        void StartOptimization() ITK_OVERRIDE;

        /** Tells Nlopt to stop the optimization at the next iteration and to returns  the best point found so far. */
        void StopOptimization()
        { nlopt_force_stop(m_NloptOptions); }

        itkGetMacro(VerboseLevel, unsigned int);
        itkSetMacro(VerboseLevel, unsigned int);

        void AddInequalityConstraint(ConstraintsFunctionType *constraint);
        void ClearInequalityConstraints();
        void AddEqualityConstraint(ConstraintsFunctionType *constraint);
        void ClearEqualityConstraints();


    protected:
        NLOPTOptimizers();
        NLOPTOptimizers(const NLOPTOptimizers&);
        virtual ~NLOPTOptimizers();
        void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;
        static double NloptFunctionWrapper(unsigned n, const double *x, double *grad, void *data);
        itkSetMacro(CurrentCost, double);

    private:
        nlopt_opt			m_NloptOptions;
        nlopt_opt           m_NloptLocalOptions;

        nlopt_algorithm		m_Algorithm;
        nlopt_algorithm     m_LocalOptimizer;

        nlopt_result		m_ErrorCode;

        ParametersType		m_LowerBoundParameters;
        ParametersType		m_UpperBoundParameters;

        bool				m_Maximize;
        bool				m_ForceStop;
        bool				m_StopValSet;
        double				m_StopVal;
        double				m_FTolRel;
        double				m_FTolAbs;
        double				m_XTolRel;
        double				m_XTolAbs;
        int                 m_MaxEval;
        double				m_MaxTime;
        int        m_VectorStorageSize;
        int        m_PopulationSize;

        unsigned int		m_VerboseLevel;

        /** Internal storage for the value type / used as a cache  */
        MeasureType			m_CurrentCost;

        std::vector<ConstraintsFunctionType::Pointer> m_InequalityConstraints;
        std::vector<ConstraintsFunctionType::Pointer> m_EqualityConstraints;

    }; // end of class

} // end of namespace anima



