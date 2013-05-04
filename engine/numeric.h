#ifndef __NUMERIC_H__
#define __NUMERIC_H__

/*! \file numeric.h
 * \brief Base file for numerical methods
 *
 * Provides basis for implementation of algorithms for numerical integration.
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

StateVector operator+(const StateVector&, const StateVector&);	///< Addition operator for StateVectors
StateVector operator-(const StateVector&, const StateVector&);	///< Subtraction operator for StateVectors
StateVector operator*(const StateVector&, double);	///< Right multiplication by double for StateVectors
StateVector operator*(double, const StateVector&);	///< Left multiplication by double for StateVectors
StateVector operator/(const StateVector&, double);	///< Division by double for StateVectors
double abs(StateVector);	///< Function returning the magnitude of a StateVector

/*! \class DiffEq
 * \brief Base class for implementing differential equations.
 *
 * Class represents a differential equation of the form dy/dx = f(x, y)
 * The (x, y) forms a state vector, which is passed to the \a derivative method, corresponding to the function f.
 */
class DiffEq
{
public:
	//! Constructor
	DiffEq();
	//! Virtual destructor
	virtual ~DiffEq();
	
	//! Pure virtual function returning the derivative based on the current state.
	/*! \param v Current state.
		\return StateVector which contains the derivatives of the components of the state.
	 */
	virtual StateVector derivative(StateVector v) = 0;
};

/*! \class Integrator
 * \brief Base class implementing a numerical integrator.
 *
 * Class represents an integrator, which calculates the state of the system based on the initial state and a time step.
 * Overloading this class allows implementation of various integrating algorithms.
 */
class Integrator
{
protected:
	double stepSize;	///< Default step size
	double initStepSize;///< Step size which was given to the constructor (for resetting purposes)
public:
	//! Constructor
	/*! \param stepSize Initial step size
	 */
	Integrator(double stepSize = 0.01);
	//! Virtual destructor
	virtual ~Integrator();
	
	//! Function returning the next value of the state.
	/*!	\param state Initial state
		\param equation Differential equation to be used for propagation
		\param step Step size. Defaults to 0.0.
		\return Next value of the state vector.
	 */
	virtual StateVector next(StateVector state, DiffEq* equation, double step = 0.0) = 0;
	
	//! Sets the default step size.
	void setStepSize(double);
	//! Resets the step size to the value passed to the constructor.
	void resetStepSize();
	
	//! Returns the default step size.
	double getStepSize();
};

#endif

