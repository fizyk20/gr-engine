#include "particle.h"

Particle::Particle(SchwManifold* _m)
	: p(0)
{
	m = _m;
}

Particle::Particle(SchwManifold* _m, Point _p, vector4 _u)
	: p(_p)
{
	m = _m;
	u = _u;
}

Particle::~Particle()
{
}

int Particle::getCoordSystem()
{
	return p.getCoordSystem();
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

void Particle::propagate(double dt)
{
	vector4 k[4][2];
	vector4 y[2], new_y[2];
	int i,j;
	
	y[0] = p.toVector4();
	y[1] = u;
	
	for(j=0; j<4; j++)
	{
		//generate k[j]
		for(i=0; i<2; i++)
		{
			new_y[i] = y[i];
			if(j>0 && j<3) new_y[i] += k[j-1][i]*dt*0.5;
			if(j==3) new_y[i] += k[j-1][i]*dt;
		}
		
		k[j][0] = new_y[1];
		Metric* metric = m->getMetric(p.getCoordSystem());
		k[j][1] = vector4() - metric->christoffel(new_y[1], new_y[1], Point(p.getCoordSystem(), new_y[0]));
	}
	
	double coeff[4] = {1.0, 2.0, 2.0, 1.0};
	for(j=0; j<4; j++)
		for(i=0; i<2; i++)
			y[i] += coeff[j]*k[j][i]*dt/6.0;
		
	p = Point(p.getCoordSystem(), y[0]);
	u = y[1];
	
	double l;
	switch(p.getCoordSystem())
	{
	case SchwManifold::EF:
		if(p[2] < 0.5)
		{
			u = m->convertVectorTo(u, p, SchwManifold::NearPole0);
			p = m->convertPointTo(p, SchwManifold::NearPole0);
			//printf("EF->Pole 0\n");
		}
		else if(p[2] > 2.642)
		{
			u = m->convertVectorTo(u, p, SchwManifold::NearPolePi);
			p = m->convertPointTo(p, SchwManifold::NearPolePi);
			//printf("EF->Pole Pi\n");
		}
		break;
	case SchwManifold::NearPole0:
	case SchwManifold::NearPolePi:
		l = p[2]*p[2] + p[3]*p[3];
		if(l > 0.07)
		{
			u = m->convertVectorTo(u, p, SchwManifold::EF);
			p = m->convertPointTo(p, SchwManifold::EF);
		}
		break;
	}
}
