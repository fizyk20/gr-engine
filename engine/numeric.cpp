#include "numeric.h"
#include <math.h>

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
		
	unsigned i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i] + arg2[i]);
		
	return result;
}

StateVector operator-(const StateVector& arg1, const StateVector& arg2)
{
	if(arg1.size() != arg2.size())
		throw StateLengthError();
		
	unsigned i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i] - arg2[i]);
		
	return result;
}

StateVector operator*(const StateVector& arg1, double arg2)
{	
	unsigned i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i]*arg2);
		
	return result;
}

StateVector operator*(double arg, const StateVector& arg2)
{	
	unsigned i;
	StateVector result;
	for(i = 0; i < arg2.size(); i++)
		result.push_back(arg2[i]*arg);
		
	return result;
}

StateVector operator/(const StateVector& arg1, double arg2)
{	
	unsigned i;
	StateVector result;
	for(i = 0; i < arg1.size(); i++)
		result.push_back(arg1[i]/arg2);
		
	return result;
}

double abs(StateVector arg)
{
	unsigned i;
	double result = 0.0;
	for(i = 0; i < arg.size(); i++)
		result += arg[i]*arg[i];
		
	return sqrt(result);
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
 
Integrator::Integrator(double stepSize)
{
	this->stepSize = stepSize;
	initStepSize = stepSize;
}

Integrator::~Integrator()
{
}

void Integrator::setStepSize(double step)
{
	stepSize = step;
}

void Integrator::resetStepSize()
{
	stepSize = initStepSize;
}

double Integrator::getStepSize()
{
	return stepSize;
}

