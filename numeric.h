#ifndef __NUMERIC_H__
#define __NUMERIC_H__

/*! \file numeric.h
 * \brief Classes implementing numerical methods
 */

#include <vector>
#include <exception>

/*! \class StateLengthError
 * \brief An exception class thrown when there is a state vector length mismatch
 */
class StateLengthError : public std::exception
{
public:
	StateLengthError();
	virtual ~StateLengthError() throw();
	const char* what() const throw();
};

/** \typedef StateVector - alias for std::vector<double> */
typedef std::vector<double> StateVector;

StateVector operator+(const StateVector&, const StateVector&);
StateVector operator-(const StateVector&, const StateVector&);
StateVector operator*(const StateVector&, double);
StateVector operator*(double, const StateVector&);
StateVector operator/(const StateVector&, double);
double abs(StateVector);

/*! \class DiffEq
 * \brief Base class for implementing differential equations.
 *        Defining differential equations is done by overloading the \a derivative method in a child class
 */
class DiffEq
{
public:
	DiffEq();
	virtual ~DiffEq();
	
	virtual StateVector derivative(StateVector v) = 0;
};

/*! \class Integrator
 * \brief Base class implementing a numerical integrator.
 */
class Integrator
{
protected:
	double stepSize;
	double initStepSize;
public:
	Integrator(double stepSize = 0.01);
	virtual ~Integrator();
	
	virtual StateVector next(StateVector state, DiffEq* equation, double step = 0.0) = 0;
	
	void setStepSize(double);
	void resetStepSize();
	
	double getStepSize();
};

#endif

