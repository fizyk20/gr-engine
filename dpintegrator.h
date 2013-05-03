#ifndef __DPINTEGRATOR__
#define __DPINTEGRATOR__

/*! \file dpintegrator.h
 * \brief Implementation of Dormand-Prince method
 */
 
 #include "numeric.h"

/*! \class DPIntegrator
 * \brief Class implementing a Dormand-Prince numerical integrator.
 */
class DPIntegrator : public Integrator
{
	double maxErr;
	double minStep;
	double maxStep;
	
	StateVector lastDerivative, lastState;
	DiffEq* lastEq;
public:
	DPIntegrator(double maxErr = 0.000001, double stepSize = 0.01);
	~DPIntegrator();
	StateVector next(StateVector state, DiffEq* equation, double step = 0.0);
	
	double getMaxErr();
	double getMinStep();
	double getMaxStep();
	
	void setMaxErr(double);
	void setMinStep(double);
	void setMaxStep(double);
};

#endif

