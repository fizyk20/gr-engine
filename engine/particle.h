#ifndef __PARTICLE__
#define __PARTICLE__

/*! \file particle.h
 * \brief Header for the Particle class, representing a point particle.
 */

#include "geometry.h"
#include "numeric.h"

/*! \class Particle
 * \brief Class representing a particle with defined position and 4-velocity.
 * 		  
 * Inherits DiffEq - defines the geodesic equation.
 */
class Particle : public DiffEq
{
protected:
	Point p;
	vector4 u;
	Manifold* m;
	
	//! Constructs a state vector from the internal state.
	virtual StateVector constructState();
	//! Sets the internal state to a state represented by a StateVector.
	virtual void setState(StateVector);
	//! Reads the position from a StateVector.
	virtual Point getPosFromState(StateVector);
	//! Reads the 4-velocity from a StateVector.
	virtual vector4 getVelFromState(StateVector);
	
	Integrator* integrator;
public:
	//! Constructor
	/*! \param _m The manifold on which the particle is defined
	 */
	Particle(Manifold* _m);
	//! Constructor
	/*! \param _m The manifold on which the particle is defined
	 *  \param _p The initial position
	 *  \param _u The initial 4-velocity
	 */
	Particle(Manifold* _m, Point _p, vector4 _u);
	//! Destructor
	virtual ~Particle();
	
	//! Sets the integrator to be used for propagation.
	void setIntegrator(Integrator*);
	
	//! Overloaded method from \a DiffEq
	/*! \param v Current state
	 *  \return The derivative of the current state as given by the geodesic equation.
	 */
	StateVector derivative(StateVector v);
	//! Propagates the particle
	/*! \param step The simulation step - corresponds to the change in proper time.
	 */
	void propagate(double step = 0.0);
	
	//! Returns the current coordinate system in use.
	int getCoordSystem();
	//! Changes the coordinate system in use.
	virtual void setCoordSystem(int);
	
	//! Returns the position.
	Point getPos();
	//! Returns the 4-velocity.
	vector4 getVel();
	
	//! Changes the position and 4-velocity
	/*! \param _p The new position
	 *  \param _u The new 4-velocity
	 */
	virtual void setPosVel(Point _p, vector4 _u);
	//! Changes the 4-velocity
	/*! \param _u The new 4-velocity
	 */
	virtual void setVel(vector4 _u);
};

#endif

