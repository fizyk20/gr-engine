#include "entity.h"
#include <math.h>

Entity::Entity(Manifold* _m, Point _p, vector4 _u, vector4 _x, vector4 _y, vector4 _z)
	: Particle(_m, _p, _u)
{
	m = _m;
	
	p = _p;
	u = _u;
	basis[0] = _x;
	basis[1] = _y;
	basis[2] = _z;
	
	orthonormalize();
}

void Entity::orthonormalize()
{
	int i,j;
	Metric* g = m->getMetric(p.getCoordSystem());
	
	vector4 v[4];
	v[0] = u;
	for(i=1; i<4; i++)
		v[i] = basis[i-1];
	
	for(i=0; i<4; i++)
	{
		for(j=0; j<i; j++)
			v[i] -= g->g(v[i], v[j], p) * v[j] / g->g(v[j], v[j], p);
		v[i] /= sqrt(fabs(g->g(v[i], v[i], p)));
	}
	
	u = v[0];
	for(i=1; i<4; i++)
		basis[i-1] = v[i];
}

Entity::~Entity()
{
}
	
StateVector Entity::constructState()
{
	StateVector result;
	int i,j;
	
	for(j=0; j<4; j++)
		result.push_back(p[j]);
	
	for(j=0; j<4; j++)
		result.push_back(u[j]);

	for(i=0; i<3; i++)
		for(j=0; j<4; j++)
			result.push_back(basis[i][j]);
			
	return result;
}

void Entity::setState(StateVector v)
{
	int i,j;
	
	for(i=0; i<4; i++)
		p[i] = v[i];
	
	for(i=0; i<4; i++)
		u[i] = v[i+4];
	
	for(j=0; j<3; j++)
		for(i=0; i<4; i++)
			basis[j][i] = v[i + 8 + j*4];
}

Point Entity::getPosFromState(StateVector v)
{
	Point result(p.getCoordSystem());
	for(int i=0; i<4; i++)
		result[i] = v[i];
		
	return result;
}

vector4 Entity::getVelFromState(StateVector v)
{
	vector4 result;
	for(int i=0; i<4; i++)
		result[i] = v[i+4];
		
	return result;
}

void Entity::setCoordSystem(int sys)
{
	if(sys != p.getCoordSystem())
	{
		u = m->convertVectorTo(u, p, sys);
		for(int i=0; i<3; i++)
			basis[i] = m->convertVectorTo(basis[i], p, sys);
		p = m->convertPointTo(p, sys);
	}
}

vector4 Entity::getStateVector(int i)
{
	if(i == 0) return u;
	else return basis[i-1];
}

vector4 Entity::getX()
{
	return basis[0];
}

vector4 Entity::getY()
{
	return basis[1];
}

vector4 Entity::getZ()
{
	return basis[2];
}

void Entity::setPosVel(Point _p, vector4 _u)
{
	p = _p;
	u = _u;
	orthonormalize();
}

void Entity::setVel(vector4 _u)
{
	u = _u;
	orthonormalize();
}

StateVector Entity::derivative(StateVector v)
{
	StateVector result;
	int i;
	for(i=0; i<4; i++)
		result.push_back(u[i]);
		
	for(i=0; i<16; i++)
		result.push_back(0.0);	//TODO
		
	return result;
}

