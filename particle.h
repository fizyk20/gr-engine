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
	
	StateVector constructState();
	void setState(StateVector);
	Point getPosFromState(StateVector);
	vector4 getVelFromState(StateVector);
	
	Integrator* integrator;
public:
	Particle(Manifold*);
	Particle(Manifold*, Point, vector4);
	virtual ~Particle();
	
	void setIntegrator(Integrator*);
	
	StateVector derivative(StateVector v);
	void propagate(double step = 0.0);
	
	int getCoordSystem();
	Point getPos();
	vector4 getVel();
	
	void setPosVel(Point, vector4);
	void setVel(vector4);
};

#endif

