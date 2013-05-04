#include "rk4integrator.h"

/*******************************************************************************
 *
 *  RK4Integrator class implementation
 *
 *******************************************************************************/

RK4Integrator::RK4Integrator(double stepSize)
	: Integrator(stepSize)
{
}

RK4Integrator::~RK4Integrator()
{
}

StateVector RK4Integrator::next(StateVector state, DiffEq* equation, double step)
{
	double h;
	if(step == 0.0) 
		h = stepSize;
	else
		h = step;
		
	StateVector k1, k2, k3, k4;
	
	k1 = equation -> derivative(state);
	k2 = equation -> derivative(state + k1*h/2);
	k3 = equation -> derivative(state + k2*h/2);
	k4 = equation -> derivative(state + k3*h);
	
	return state + (k1 + 2*k2 + 2*k3 + k4)*h/6;
}

