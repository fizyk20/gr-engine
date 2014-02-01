#ifndef __ENTITY__
#define __ENTITY__

/*! \file entity.h
 * \brief Header for the Entity class, representing an entity having an orientation.
 */

#include "particle.h"
#include "numeric.h"

/*! \class Entity
 * \brief Class representing an entity with defined position and local basis, consisting of 4-velocity and 3 spatial directions.
 * 
 * The entity can be accelerated and rotated.
 * Inherits Particle.
 */
class Entity : public Particle
{
protected:
	vector4 basis[3]; //local basis
	Manifold* m;
	
	//! Orthonormalizes the local basis, preserving the direction of the 4-velocity (other directions may be changed).
	void orthonormalize();
	
	double force[3];
	double angvel[3];
	
	//! Calculates the covariant derivative of a component of the local basis.
	/*! \param component The number of the component (0 - 4-velocity etc.)  (\a getStateVector)
	 *	\return The covariant derivative of the selected component
	 */
	vector4 calculateFourForce(int component);
	
	//! Constructs a state vector from the internal state.
	StateVector constructState();
	//! Sets the internal state to a state represented by a StateVector.
	void setState(StateVector);
	//! Reads the position from a StateVector.
	Point getPosFromState(StateVector);
	//! Reads the 4-velocity from a StateVector.
	vector4 getVelFromState(StateVector);
	//! Reads a component of the local basis from a StateVector.
	vector4 getVectorFromState(StateVector, int);
public:
	//Entity(Manifold*, Point, vector4); - TODO
	//! Constructor
	/*! Initializes the entity with the manifold, position and local basis. The local basis is then orthonormalized (preserving the direction of the 4-velocity).
	 * \param _m The manifold to be used
	 * \param _p Initial position
	 * \param _u Initial 4-velocity
	 * \param _x Initial local X direction
	 * \param _y Initial local Y direction
	 * \param _z Initial local Z direction
	 */
	Entity(Manifold* _m, Point _p, vector4 _u, vector4 _x, vector4 _y, vector4 _z);
	//! Destructor
	~Entity();
	
	//! Returns a component of the local basis.
	/*! \param i The number of the component (0 - 4-velocity etc, 1 - local X, 2 - local Y, 3 - local Z)
	 */
	vector4 getStateVector(int i);
	//! Returns the local X direction
	vector4 getX();
	//! Returns the local Y direction
	vector4 getY();
	//! Returns the local Z direction
	vector4 getZ();
	
	//! Overloaded method from \a DiffEq
	/*! \param v Current state
	 *  \return The derivative of the current state as given by the Fermi-Walker transport.
	 */
	StateVector derivative(StateVector v);
	//! Overloaded method changing the coordinate system in use.
	void setCoordSystem(int);
	
	//! Changes the position and 4-velocity
	/*! \param _p The new position
	 *  \param _u The new 4-velocity
	 */
	void setPosVel(Point, vector4);
	//! Changes the 4-velocity
	/*! \param _u The new 4-velocity
	 */
	void setVel(vector4);
	
	//! Adds an acceleration to the entity
	/*! \param x Local X component
	 *  \param y Local Y component
	 *  \param z Local Z component
	 */
	void applyForce(double x, double y, double z);
	//! Adds angular velocity to the entity
	/*! \param x Local X component
	 *  \param y Local Y component
	 *  \param z Local Z component
	 */
	void applyAngVel(double x, double y, double z);
	//! Rotates the entity
	/*! \param pitch Rotation pitch angle
	 *  \param yaw Rotation yaw angle
	 *  \param roll Rotation roll angle
	 */
	void rotate(double pitch, double yaw, double roll);
	
	//! Overloaded Particle::propagate, zeroing force and ang-vel after each step
	void propagate(double step=0.0);
};

#endif

