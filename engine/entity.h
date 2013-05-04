#ifndef __ENTITY__
#define __ENTITY__

/*! \file entity.h
 * \brief Header for the Entity class, representing an entity having an orientation.
 */

#include "particle.h"
#include "numeric.h"

/*! \class Entity
 * \brief Class representing an entity with defined position, 4-velocity and local basis (orientation)
 * 		  Inherits Particle
 */
class Entity : public Particle
{
	Point p;
	vector4 basis[3]; //local basis
	Manifold* m;
	
	void orthonormalize();
	
	double force[3];
	double angvel[3];
	
	vector4 calculateFourForce(int);
	
	StateVector constructState();
	void setState(StateVector);
	Point getPosFromState(StateVector);
	vector4 getVelFromState(StateVector);
	vector4 getVectorFromState(StateVector, int);
public:
	//Entity(Manifold*, Point, vector4); - TODO
	Entity(Manifold*, Point, vector4, vector4, vector4, vector4);
	~Entity();
	
	vector4 getStateVector(int);
	vector4 getX();
	vector4 getY();
	vector4 getZ();
	
	StateVector derivative(StateVector v);
	void setCoordSystem(int);
	
	void setPosVel(Point, vector4);
	void setVel(vector4);
	
	void applyForce(double, double, double);
	void applyAngVel(double, double, double);
	void rotate(double, double, double);
};

#endif
