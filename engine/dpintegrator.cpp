#include "dpintegrator.h"
#include <math.h>

/*******************************************************************************
 *
 *  DPIntegrator class implementation
 *
 *******************************************************************************/

DPIntegrator::DPIntegrator(double maxErr, double stepSize, double minStep, double maxStep)
	: Integrator(stepSize)
{
	this->maxErr = maxErr;
	this->minStep = minStep;
	this->maxStep = maxStep;
}

DPIntegrator::~DPIntegrator()
{
}

StateVector DPIntegrator::next(StateVector state, DiffEq* equation, double step)
{
	double h;
	if(step == 0.0) 
		h = stepSize;
	else
		h = step;
	
	StateVector k1, k2, k3, k4, k5, k6, k7, nextState;
	
	//optimization - use only if last calculated derivative exists, if the same equation is used and if neither state nor time were changed
	if(lastDerivative.size() && lastEq == equation && lastState == state)
		k1 = h * lastDerivative;
	else
		k1 = h * equation -> derivative(state);
	k2 = h * equation -> derivative(state + k1/5);
	k3 = h * equation -> derivative(state + k1*3.0/40 + k2*9.0/40);
	k4 = h * equation -> derivative(state + k1*44.0/45 - k2*56.0/15 + k3*32.0/9);
	k5 = h * equation -> derivative(state + k1*19372.0/6561 - k2*25360.0/2187 + k3*64448.0/6561 - k4*212.0/729);
	k6 = h * equation -> derivative(state + k1*9017.0/3168 - k2*355.0/33 + k3*46732.0/5247 + k4*49.0/176 - k5*5103.0/18656);


	nextState = state + k1*35.0/384 + k3*500.0/1113 + k4*125.0/192 - k5*2187.0/6784 + k6*11.0/84;

	k7 = equation -> derivative(nextState);

	double error = abs(k1*71.0/57600 - k3*71.0/16695 + k4*71.0/1920 - k5*17253.0/339200 + k6*22.0/525 - h*k7/40);
	
	if(error != 0.0) stepSize = h*pow(maxErr/error, 0.25);
	else stepSize = maxStep;
	
	if(stepSize < minStep) stepSize = minStep;
	if(stepSize > maxStep) stepSize = maxStep;
	if(stepSize < 0.8*h && step == 0.0)
		return next(state, equation);
	
	//for optimization
	lastDerivative = k7;
	lastState = nextState;
	lastEq = equation;
		
	return nextState;
}

double DPIntegrator::getMaxErr()
{
	return maxErr;
}

double DPIntegrator::getMinStep()
{
	return minStep;
}

double DPIntegrator::getMaxStep()
{
	return maxStep;
}

void DPIntegrator::setMaxErr(double mE)
{
	maxErr = mE;
}

void DPIntegrator::setMinStep(double mS)
{
	minStep = mS;
}

void DPIntegrator::setMaxStep(double mS)
{
	maxStep = mS;
}

