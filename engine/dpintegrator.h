#ifndef __DPINTEGRATOR__
#define __DPINTEGRATOR__

/*! \file dpintegrator.h
 * \brief Implementation of Dormand-Prince method
 */
 
 #include "numeric.h"

/*! \class DPIntegrator
 * \brief Class implementing a Dormand-Prince numerical integrator.
 *
 * The Dormand-Prince method is an adaptive method. After each step the step size is recalculated, so that the error will be just below the predefined error margin. The minimal and maximal step size can also be set.
 */
class DPIntegrator : public Integrator
{
	double maxErr;
	double minStep;
	double maxStep;
	
	StateVector lastDerivative, lastState;
	DiffEq* lastEq;
public:
	//! Constructor
	/*! \param maxErr The error margin - if the error is larger than this margin, the step size is decreased.
	 *  \param stepSize Default step size
	 *  \param minStep Minimal step size
	 *  \param maxStep Maximal step size
	 */
	DPIntegrator(double maxErr = 0.000001, double stepSize = 0.01, double minStep = 0.0001, double maxStep = 0.1);
	//! Destructor
	~DPIntegrator();
	//! Function calculating the next state
	/*! \param state Current state
	 *  \param equation The differential equation to be used
	 *  \param step Step size. If 0 (default), the default step size is used.
	 */
	StateVector next(StateVector state, DiffEq* equation, double step = 0.0);
	
	//! Returns the error margin
	double getMaxErr();
	//! Returns the minimal step size
	double getMinStep();
	//! Returns the maximal step size
	double getMaxStep();
	
	//! Sets the error margin
	void setMaxErr(double);
	//! Sets the minimal step size
	void setMinStep(double);
	//! Sets the maximal step size
	void setMaxStep(double);
};

#endif

