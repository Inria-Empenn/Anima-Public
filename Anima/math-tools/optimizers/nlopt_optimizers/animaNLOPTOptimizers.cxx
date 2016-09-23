/*
 * ITK Wrapper for the NLOPT library
 *
 * This code was graciously provided to us by the CRL at Boston. Please do NOT distribute
 * it in any form or to anyone.
 *
 * Copyright (c) 2010-2011 Children's Hospital Boston.
 * Benoit Scherrer, CRL (Computational Radiology Laboratory), Harvard Medical School
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/
#include "animaNLOPTOptimizers.h"

#include <math.h>
#include <vnl/vnl_math.h>

#include <itkTimeProbe.h>

using namespace std;

namespace anima
{
    /**********************************************************************************************//**
     * \fn	NLOPTOptimizers::NLOPTOptimizers()
     *
     * \brief	Default constructor.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
    *************************************************************************************************/
    NLOPTOptimizers::NLOPTOptimizers()
    {
        m_VerboseLevel = 0;
        m_Maximize = false;
        m_ErrorCode = NLOPT_FAILURE;
        m_Algorithm = NLOPT_LN_BOBYQA;
        m_LocalOptimizer = NLOPT_LN_BOBYQA;

        m_StopValSet = false;
        m_StopVal = 0;
        m_FTolRel = -1;
        m_FTolAbs = -1;
        m_XTolRel = -1;
        m_XTolAbs = -1;
        m_MaxEval = -1;
        m_MaxTime = -1;
        m_PopulationSize = -1;
        m_VectorStorageSize = -1;
        m_CurrentCost = 0;
    }

    /**********************************************************************************************//**
     * \fn	NLOPTOptimizers::~NLOPTOptimizers()
     *
     * \brief	Destructor.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
    *************************************************************************************************/
    NLOPTOptimizers::~NLOPTOptimizers()
    {
    }

    /**********************************************************************************************//**
     * \fn	double NLOPTOptimizers::NloptFunctionWrapper(unsigned n, const double *x, double *grad, void *data)
     *
     * \brief	NloptFunctionWrapper is called by nlopt to evaluate the function f to minimize.
     * 			NloptFunctionWrapper converts the nlopt parameters to the SingleValuedCostFunction
     * 			parameters, and call the SingleValuedCostFunction to actually compute the value of f.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \param	n				The.
     * \param	x				The x coordinate.
     * \param [in,out]	grad	If non-null, the graduated.
     * \param [in,out]	data	If non-null, the data.
     *
     * \return	.
    *************************************************************************************************/
    double NLOPTOptimizers::NloptFunctionWrapper(unsigned n, const double *x, double *grad, void *data)
    {
        NLOPTOptimizers *optimizer = static_cast<NLOPTOptimizers *>(data);

        //-----------------------------------------
        // Copy the nlopt position to a itk position
        // (takes into account the itk scale)
        //-----------------------------------------
        NLOPTOptimizers::ParametersType itkCurrentPosition(n);
        for ( unsigned int i=0; i<n ; i++ )
            itkCurrentPosition[i] = x[i]/optimizer->GetScales()[i];

        double f = (optimizer->GetCostFunction()->GetValue(itkCurrentPosition));

        //-----------------------------------------
        // If needed, also compute the derivative and copy
        // back in *grad
        //-----------------------------------------
        if ( grad!=NULL )
        {
            DerivativeType derivative;
            optimizer->GetCostFunction()->GetDerivative (itkCurrentPosition, derivative);
            for ( unsigned int i=0; i<n ; i++ )
                grad[i] = derivative[i];
        }

        return f;
    }

    /**********************************************************************************************//**
     * \fn	void NLOPTOptimizers::StartOptimization()
     *
     * \brief	Starts the optimization.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \exception	itk::ExceptionObject	Thrown when exception.
    *************************************************************************************************/
    void NLOPTOptimizers::StartOptimization()
    {
        if( m_CostFunction.IsNull() )
            throw itk::ExceptionObject(__FILE__, __LINE__, "Error. The cost function has not been set", "NLOPTOptimizers");

        //---------------------------------------------
        // n = number of parameters
        //---------------------------------------------
        unsigned int n = m_CostFunction->GetNumberOfParameters();

        //---------------------------------------------
        // Set x = initial position
        // Take into account the ITK scales
        //---------------------------------------------
        if ( this->GetInitialPosition().GetSize() != n )
                throw itk::ExceptionObject(__FILE__, __LINE__, "Invalid initial position parameter. Its size should be equal to the number of parameters", "NLOPTOptimizers");
        m_CurrentPosition.set_size(n);

        double *x = new double[n];
        const double *in_x = GetInitialPosition().data_block();
        for ( unsigned int i=0; i<n ; i++ )
        {
            x[i] = in_x[i] * GetScales()[i];
            m_CurrentPosition[i] = in_x[i];
        }

        //---------------------------------------------
        // lb/ub = lower bound/upper bound
        // Take into account the ITK scales
        //---------------------------------------------
        double *lb=NULL;
        double *ub=NULL;
        if ( m_LowerBoundParameters.GetSize()!=0 )
        {
            if ( m_LowerBoundParameters.GetSize()!=n )
                throw itk::ExceptionObject(__FILE__, __LINE__, "Invalid lower bound parameter. Its size should be equal to the number of parameters", "NLOPTOptimizers");

            lb = new double[n];
            for ( unsigned int i=0; i<n; i++ )
                lb[i] = m_LowerBoundParameters[i]*GetScales()[i];
        }
        if ( m_UpperBoundParameters.GetSize()!=0 )
        {
            if ( m_UpperBoundParameters.GetSize()!=n )
                throw itk::ExceptionObject(__FILE__, __LINE__, "Invalid upper bound parameter. Its size should be equal to the number of parameters", "NLOPTOptimizers");

            ub = new double[n];
            for ( unsigned int i=0; i<n; i++ )
                ub[i] = m_UpperBoundParameters[i]*GetScales()[i];
        }
        
        //---------------------------------------------
        // Creates the NLOPT structure and fill it
        //---------------------------------------------
        m_NloptOptions = nlopt_create((::nlopt_algorithm)(int)m_Algorithm, n);
        nlopt_set_lower_bounds(m_NloptOptions, lb);
        nlopt_set_upper_bounds(m_NloptOptions, ub);
        if ( m_Maximize )
            nlopt_set_max_objective(m_NloptOptions, (nlopt_func)this->NloptFunctionWrapper, (void *)this);
        else
            nlopt_set_min_objective(m_NloptOptions, (nlopt_func)this->NloptFunctionWrapper, (void *)this);
        
        if ( m_StopValSet ) nlopt_set_stopval(m_NloptOptions, m_StopVal);
        
        nlopt_set_ftol_rel(m_NloptOptions, m_FTolRel);
        nlopt_set_ftol_abs(m_NloptOptions, m_FTolAbs);
        nlopt_set_xtol_rel(m_NloptOptions, m_XTolRel);
        nlopt_set_xtol_abs1(m_NloptOptions, m_XTolAbs);
        nlopt_set_maxtime(m_NloptOptions, m_MaxTime);
        nlopt_set_vector_storage(m_NloptOptions, m_VectorStorageSize);
        nlopt_set_maxeval(m_NloptOptions, m_MaxEval);
        nlopt_set_population(m_NloptOptions, m_PopulationSize);
        
        for (unsigned int i = 0;i < m_InequalityConstraints.size();++i)
            nlopt_add_inequality_constraint(m_NloptOptions, (nlopt_func)ConstraintsFunctionType::GetConstraintValue, m_InequalityConstraints[i]->GetAdditionalData(), m_InequalityConstraints[i]->GetTolerance());
        
        for (unsigned int i = 0;i < m_EqualityConstraints.size();++i)
            nlopt_add_equality_constraint(m_NloptOptions, (nlopt_func)ConstraintsFunctionType::GetConstraintValue, m_EqualityConstraints[i]->GetAdditionalData(), m_EqualityConstraints[i]->GetTolerance());
        
        //----------------------------------------
        // Setup local optimizer for algorithms
        // using sequences of local optimizations
        //----------------------------------------
        m_NloptLocalOptions = nlopt_create((::nlopt_algorithm)(int)m_LocalOptimizer, n);
        
        if ( m_StopValSet ) nlopt_set_stopval(m_NloptLocalOptions, m_StopVal);
        
        nlopt_set_ftol_rel(m_NloptLocalOptions, m_FTolRel);
        nlopt_set_ftol_abs(m_NloptLocalOptions, m_FTolAbs);
        nlopt_set_xtol_rel(m_NloptLocalOptions, m_XTolRel);
        nlopt_set_xtol_abs1(m_NloptLocalOptions, m_XTolAbs);
        nlopt_set_maxtime(m_NloptLocalOptions, m_MaxTime);
        nlopt_set_vector_storage(m_NloptLocalOptions, m_VectorStorageSize);
        nlopt_set_maxeval(m_NloptLocalOptions, m_MaxEval);
        
        //--------------------------------------
        // Plug local optimizer into global one
        //--------------------------------------
        nlopt_set_local_optimizer(m_NloptOptions, m_NloptLocalOptions);

        this->InvokeEvent( itk::StartEvent() );

        //----------------------------------------
        // Run the NLOPT optimizer !
        //----------------------------------------
        double valf;
        m_ErrorCode = static_cast<nlopt_result>((int)nlopt_optimize (m_NloptOptions, x, &valf));
        SetCurrentCost(valf);

        //----------------------------------------
        // Converts back using the scales
        //----------------------------------------
        NLOPTOptimizers::ParametersType p(n);
        for ( unsigned int i=0; i<n ; i++ )
            p[i] = x[i]/GetScales()[i];
        this->SetCurrentPosition(p);

        //----------------------------------------
        // Free memory
        //----------------------------------------
        delete[] x;
        if ( ub!=NULL ) delete[] ub;
        if ( lb!=NULL ) delete[] lb;

        nlopt_destroy(m_NloptOptions);
        nlopt_destroy(m_NloptLocalOptions);

        this->InvokeEvent( itk::EndEvent() );
    }

    /**********************************************************************************************//**
     * \fn	bool NLOPTOptimizers::isSuccessful() const
     *
     * \brief	Return true if the optimization was successful (ie the last error code was positive
     * 			or equal to NLOPT_ROUNDOFF_LIMITED or NLOPT_FORCED_STOP)
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \return	true if successful, false if not.
    *************************************************************************************************/
    bool NLOPTOptimizers::isSuccessful() const
    {
        if (m_ErrorCode < 0
            && m_ErrorCode != NLOPT_ROUNDOFF_LIMITED
            && m_ErrorCode != NLOPT_FORCED_STOP)
            return false;
        else
            return true;
    }

    /**********************************************************************************************//**
     * \fn	std::string NLOPTOptimizers::GetErrorCodeDescription() const
     *
     * \brief	Gets the current error code description as a std::string.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \return	The error code description.
    *************************************************************************************************/
    std::string NLOPTOptimizers::GetErrorCodeDescription() const
    {
        switch ( m_ErrorCode )
        {
        case NLOPT_SUCCESS:
            return "Success";

        case NLOPT_STOPVAL_REACHED:
            return "Optimization stopped because 'StopVal' was reached.";
        case NLOPT_FTOL_REACHED:
            return "Optimization stopped because 'FTolRel' or 'FTolAbs' was reached.";
        case NLOPT_XTOL_REACHED:
            return "Optimization stopped because 'XTolRel' or 'XTolAbs' was reached.";
        case NLOPT_MAXEVAL_REACHED:
            return "Optimization stopped because 'MaxEval' was reached.";
        case NLOPT_MAXTIME_REACHED:
            return "Optimization stopped because 'MaxTime' was reached.";


        case NLOPT_FAILURE:
            return "Generic failure code.";
        case NLOPT_INVALID_ARGS:
            return "Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).";
        case NLOPT_OUT_OF_MEMORY:
            return "Out of memory.";
        case NLOPT_ROUNDOFF_LIMITED:
            return "Halted because roundoff errors limited progress.";
        case NLOPT_FORCED_STOP:
            return "Halted because of a forced termination.";

        default:
            break;
        }
        return "Unkown error code";
    }

    /**********************************************************************************************//**
     * \fn	void NLOPTOptimizers::AddInequalityConstraint(ConstraintsFunctionType::Pointer &constraint)
     *
     * \brief	Adds a nonlinear inequality constraint of the form \f$ fc(x) <= 0 \f$. The parameter
     * 			tol is a tolerance that is used for the purpose of stopping criteria only: a point x
     * 			is considered feasible for judging whether to stop the optimization if \f$ fc(x) <=
     * 			tol \f$. A tolerance of zero means that NLopt will try not to consider any x to be
     * 			converged unless fc is strictly non-positive; generally, at least a small positive
     * 			tolerance is advisable to reduce sensitivity to rounding errors.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \param	constraint The constraint function containing all details
    *************************************************************************************************/
    void NLOPTOptimizers::AddInequalityConstraint(ConstraintsFunctionType *constraint)
    {
        m_InequalityConstraints.push_back(constraint);
    }

    /**********************************************************************************************//**
     * \fn	void NLOPTOptimizers::ClearInequalityConstraints()
     *
     * \brief	Clears all inequality constraints.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
    *************************************************************************************************/
    void NLOPTOptimizers::ClearInequalityConstraints()
    {
        m_InequalityConstraints.clear();
    }

    /**********************************************************************************************//**
     * \fn	void NLOPTOptimizers::AddEqualityConstraint(ConstraintsFunctionType::Pointer &constraint)
     *
     * \brief	Adds a nonlinear equality constraint of the form \f$ h(x) = 0 \f$. The parameter tol
     * 			is a tolerance that is used for the purpose of stopping criteria only: a point x is
     * 			considered feasible for judging whether to stop the optimization if \f$ |h(x)| <= tol
     * 			\f$. For equality constraints, a small positive tolerance is strongly advised in
     * 			order to allow NLopt to converge even if the equality constraint is slightly nonzero.
     *
     * 			(For any algorithm listed as "derivative-free" below, the grad argument to fc or h
     * 			will always be NULL and need never be computed.)
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \param	constraint The constraint function containing all details
    *************************************************************************************************/
    void NLOPTOptimizers::AddEqualityConstraint(ConstraintsFunctionType *constraint)
    {
        m_EqualityConstraints.push_back(constraint);
    }

    /**********************************************************************************************//**
     * \fn	void NLOPTOptimizers::ClearEqualityConstraints()
     *
     * \brief	Clears all equality constraints.
     *
     * \author	Benoit Scherrer
     * \date	July 2010
    *************************************************************************************************/
    void NLOPTOptimizers::ClearEqualityConstraints()
    {
        m_EqualityConstraints.clear();
    }

    /**********************************************************************************************//**
     * \fn	void NLOPTOptimizers::SetLowerBoundParameters( const ParametersType& p )
     *
     * \brief	Set the lower bound parameters (not used by all algorithms)
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \param	p	The lower bound parameters.
    *************************************************************************************************/
    void NLOPTOptimizers::SetLowerBoundParameters( const ParametersType& p )
    {
        m_LowerBoundParameters.SetSize(p.GetSize());
        for ( unsigned int i=0; i<p.GetSize(); i++ )
            m_LowerBoundParameters[i] = p[i];
    }

    /**********************************************************************************************//**
     * \fn	void NLOPTOptimizers::SetUpperBoundParameters( const ParametersType& p )
     *
     * \brief	Set the upper bound parameters (not used by all algorithms)
     *
     * \author	Benoit Scherrer
     * \date	July 2010
     *
     * \param	p	The upper bound parameters.
    *************************************************************************************************/
    void NLOPTOptimizers::SetUpperBoundParameters( const ParametersType& p )
    {
        m_UpperBoundParameters.SetSize(p.GetSize());
        for ( unsigned int i=0; i<p.GetSize(); i++ )
            m_UpperBoundParameters[i] = p[i];
    }

    /**
    *
    */
    void NLOPTOptimizers::PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);

        int major, minor, bugfix;
        nlopt_version(&major, &minor, &bugfix);
        os << indent << "NLOPT v"<<major<<"."<<minor<<"."<<bugfix<<std::endl;

        os << indent << "Algorithm: " << ::nlopt_algorithm_name((::nlopt_algorithm)(int)m_Algorithm) << std::endl;
        os << indent << "Maximize On/Off   " << m_Maximize << std::endl;

        os << indent << "StopVal: "<<m_StopVal<<std::endl;
        os << indent << "FTolRel: "<<m_FTolRel<<std::endl;
        os << indent << "FTolAbs: "<<m_FTolAbs<<std::endl;
        os << indent << "XTolRel: "<<m_XTolRel<<std::endl;
        os << indent << "MaxEval: "<<m_MaxEval<<std::endl;
        os << indent << "MaxTime: "<<m_MaxTime<<std::endl;
    }

} // end of namespace anima
