#ifndef __RK4INTEGRATOR__
#define __RK4INTEGRATOR__

/*! \file rk4integrator.h
 * \brief Implementation of RK4 method
 */
 
 #include "numeric.h"

/*! \class RK4Integrator
 * \brief Class implementing an RK4 numerical integrator.
 */
class RK4Integrator : public Integrator
{
public:
	RK4Integrator(double stepSize = 0.01);
	~RK4Integrator();
	StateVector next(StateVector state, DiffEq* equation, double step = 0.0);
};

#endif

