#include "particle.h"

Particle::Particle(Manifold* _m)
	: p(0)
{
	m = _m;
	integrator = NULL;
}

Particle::Particle(Manifold* _m, Point _p, vector4 _u)
	: p(_p)
{
	m = _m;
	u = _u;
	integrator = NULL;
}

Particle::~Particle()
{
	integrator = NULL;
}

int Particle::getCoordSystem()
{
	return p.getCoordSystem();
}

void Particle::setCoordSystem(int sys)
{
	if(sys != p.getCoordSystem())
	{
		u = m->convertVectorTo(u, p, sys);
		p = m->convertPointTo(p, sys);
	}
}

Point Particle::getPos()
{
	return p;
}

vector4 Particle::getVel()
{
	return u;
}

void Particle::setPosVel(Point _p, vector4 _u)
{
	p = _p;
	u = _u;
}

void Particle::setVel(vector4 _u)
{
	u = _u;
}

StateVector Particle::constructState()
{
	StateVector v;
	int i;
	for(i = 0; i < 4; i++)
	{
		v.push_back(p[i]);
		v.push_back(u[i]);
	}
	return v;
}

void Particle::setState(StateVector v)
{
	int i, j;
	j = 0;
	for(i = 0; i < 4; i++)
	{
		p[i] = v[j];
		j++;
		u[i] = v[j];
		j++;
	}
}

Point Particle::getPosFromState(StateVector v)
{
	Point p1(p.getCoordSystem());
	int i;
	for(i = 0; i < 4; i++)
		p1[i] = v[2*i];
		
	return p1;
}

vector4 Particle::getVelFromState(StateVector v)
{
	vector4 v1;
	int i;
	for(i = 0; i < 4; i++)
		v1[i] = v[2*i+1];
		
	return v1;
}

void Particle::setIntegrator(Integrator* i)
{
	integrator = i;
}

void Particle::propagate(double dt)
{
	if(!integrator) throw "Integrator not set!";
	
	setState(integrator -> next(constructState(), this, dt));
	
	int newCoordSystem = m->recommendCoordSystem(p);
	setCoordSystem(newCoordSystem);
}

StateVector Particle::derivative(StateVector v)
{
	Point p1 = getPosFromState(v);
	vector4 u1 = getVelFromState(v);
	
	Metric* metric = m -> getMetric(p.getCoordSystem());
	vector4 du = vector4() - metric -> christoffel(u1, u1, p1);
	
	StateVector result;
	int i;
	for(i = 0; i < 4; i++)
	{
		result.push_back(u1[i]);
		result.push_back(du[i]);
	}
	return result;
}

