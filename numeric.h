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

/*! \class DiffEq
 * \brief Base class for implementing differential equations.
 *        Defining differential equations is done by overloading the \a derivative method in a child class
 */
class DiffEq
{
public:
	DiffEq();
	virtual ~DiffEq();
	
	virtual StateVector derivative(double t, StateVector v) = 0;
};

/*! \class Integrator
 * \brief Base class implementing a numerical integrator.
 */
class Integrator
{
protected:
	double t;
	double stepSize;
	DiffEq* equation;
	StateVector currentState;
public:
	Integrator(double t0 = 0.0, double stepSize = 0.01);
	virtual ~Integrator();
	
	virtual StateVector next(double step = 0.0) = 0;
	
	void setStepSize(double);
	void setDiffEquation(DiffEq*);
	void setState(StateVector);
	void setT(double);
	
	double getStepSize();
	StateVector getState();
	double getT();
};

/*! \class RK4Integrator
 * \brief Class implementing an RK4 numerical integrator.
 */
class RK4Integrator : public Integrator
{
public:
	RK4Integrator(double t0 = 0.0, double stepSize = 0.01);
	~RK4Integrator();
	StateVector next(double step = 0.0);
};

#endif

