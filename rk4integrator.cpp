#include "rk4integrator.h"

/*******************************************************************************
 *
 *  RK4Integrator class implementation
 *
 *******************************************************************************/

RK4Integrator::RK4Integrator(double stepSize, double t0)
	: Integrator(stepSize, t0)
{
}

RK4Integrator::~RK4Integrator()
{
}

StateVector RK4Integrator::next(DiffEq* equation, double step)
{
	double h;
	if(step == 0.0) 
		h = stepSize;
	else
		h = step;
		
	StateVector k1, k2, k3, k4;
	
	k1 = equation -> derivative(t, currentState);
	k2 = equation -> derivative(t + h/2, currentState + k1*h/2);
	k3 = equation -> derivative(t + h/2, currentState + k2*h/2);
	k4 = equation -> derivative(t + h, currentState + k3*h);
	
	currentState = currentState + (k1 + 2*k2 + 2*k3 + k4)*h/6;
	t += h;
		
	return currentState;
}

