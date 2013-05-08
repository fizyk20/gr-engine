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
	//! Constructor
	/*! \param stepSize Default step size
	 */
	RK4Integrator(double stepSize = 0.01);
	//! Destructor
	~RK4Integrator();
	//! Function calculating the next state
	/*! \param state Current state
	 *  \param equation The differential equation to be used
	 *  \param step Step size. If 0 (default), the default step size is used.
	 */
	StateVector next(StateVector state, DiffEq* equation, double step = 0.0);
};

#endif

