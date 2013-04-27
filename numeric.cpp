#include "numeric.h"

StateLengthError::StateLengthError()
{
}

StateLengthError::~StateLengthError() throw()
{
}

// StateLengthError "what" method

const char* StateLengthError::what() const throw()
{
	return "State vector length mismatch.";
}

/*******************************************************************************
 *
 *  StateVector operators implementation
 *
 *******************************************************************************/
 
StateVector operator+(const StateVector& arg1, const StateVector& arg2)
{
	if(arg1.size() != arg2.size())
		throw StateLengthError();
		
	int i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i] + arg2[i]);
		
	return result;
}

StateVector operator-(const StateVector& arg1, const StateVector& arg2)
{
	if(arg1.size() != arg2.size())
		throw StateLengthError();
		
	int i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i] - arg2[i]);
		
	return result;
}

StateVector operator*(const StateVector& arg1, double arg2)
{	
	int i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i]*arg2);
		
	return result;
}

StateVector operator*(double arg, const StateVector& arg2)
{	
	int i;
	StateVector result;
	for(i = 0; i < arg2.size(); i++)
		result.push_back(arg2[i]*arg);
		
	return result;
}

StateVector operator/(const StateVector& arg1, double arg2)
{	
	int i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i]/arg2);
		
	return result;
}

/*******************************************************************************
 *
 *  DiffEq class implementation
 *
 *******************************************************************************/
 
DiffEq::DiffEq()
{
}

DiffEq::~DiffEq()
{
}

/*******************************************************************************
 *
 *  Integrator class implementation
 *
 *******************************************************************************/
 
Integrator::Integrator(double t0, double stepSize)
{
	t = t0;
	this->stepSize = stepSize;
}

Integrator::~Integrator()
{
}

void Integrator::setStepSize(double step)
{
	stepSize = step;
}

void Integrator::setDiffEquation(DiffEq* eq)
{
	equation = eq;
}

void Integrator::setState(StateVector v)
{
	currentState = v;
}

void Integrator::setT(double t)
{
	this->t = t;
}

double Integrator::getStepSize()
{
	return stepSize;
}

StateVector Integrator::getState()
{
	return currentState;
}

double Integrator::getT()
{
	return t;
}

/*******************************************************************************
 *
 *  RK4Integrator class implementation
 *
 *******************************************************************************/

RK4Integrator::RK4Integrator(double t0, double stepSize)
	: Integrator(t0, stepSize)
{
}

RK4Integrator::~RK4Integrator()
{
}

StateVector RK4Integrator::next(double step)
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

