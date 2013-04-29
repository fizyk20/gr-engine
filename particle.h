#ifndef __PARTICLE__
#define __PARTICLE__

/*! \file particle.h
 * \brief Header for the Particle class, representing a point particle.
 */

#include "geometry.h"
#include "numeric.h"

/*! \class Particle
 * \brief Class representing a particle with defined position and 4-velocity
 * 		  Inherits DiffEq for solving of the geodesic equation
 */
class Particle : public DiffEq
{
protected:
	Point p;
	vector4 u;
	Manifold* m;
	
	virtual StateVector constructState();
	virtual void setState(StateVector);
public:
	Particle(Manifold*);
	Particle(Manifold*, Point, vector4);
	virtual ~Particle();
	
	void propagate(double step = 0.0);
	
	int getCoordSystem();
	Point getPos();
	vector4 getVel();
	
	void setPosVel(Point, vector4);
	void setVel(vector4);
};

#endif

