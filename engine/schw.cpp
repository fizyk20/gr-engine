#include "schw.h"
#include <math.h>

SchwManifold::SchwManifold(double _M)
{
	M = _M;
	
	gManifold = this;
	Point::setGlobalManifold(this);
	
	nCoordSystems = 3;
	
	metrics = new Metric*[nCoordSystems];
	metrics[EF] = new SchwEFMetric(EF, this);
	metrics[NearPole0] = new SchwNearPoleMetric(NearPole0, this);
	metrics[NearPolePi] = new SchwNearPoleMetric(NearPolePi, this);
	
	conversions = new CoordinateConversion**[nCoordSystems];
	conversions[0] = new CoordinateConversion*[nCoordSystems];
	conversions[1] = new CoordinateConversion*[nCoordSystems];
	conversions[2] = new CoordinateConversion*[nCoordSystems];
	
	conversions[EF][EF] = new IdentityConversion;
	conversions[NearPole0][NearPole0] = new IdentityConversion;
	conversions[NearPolePi][NearPolePi] = new IdentityConversion;
	conversions[EF][NearPole0] = new EFToNearPole0;
	conversions[NearPole0][EF] = new NearPole0ToEF;
	conversions[EF][NearPolePi] = new EFToNearPolePi;
	conversions[NearPolePi][EF] = new NearPolePiToEF;
	conversions[NearPole0][NearPolePi] = new NearPoleToNearPole;
	conversions[NearPolePi][NearPole0] = new NearPoleToNearPole;
}

SchwManifold::~SchwManifold()
{
	int i,j;
	for(i=0; i<nCoordSystems; i++)
		delete metrics[i];
	delete[] metrics;
	
	for(i=0; i<nCoordSystems; i++)
	{
		for(j=0; j<nCoordSystems; j++)
			delete conversions[i][j];
		delete[] conversions[i];
	}
	delete[] conversions;
	
	nCoordSystems = 0;
}

double SchwManifold::getMass()
{
	return M;
}

void SchwManifold::setMass(double _M)
{
	M = _M;
}

int SchwManifold::recommendCoordSystem(Point p)
{
	double l;
	switch(p.getCoordSystem())
	{
	case EF:
		if(p[2] < 0.5) return NearPole0;
		else if(p[2] > 2.642) return NearPolePi;
		else return EF;
	case NearPole0:
	case NearPolePi:
		l = p[2]*p[2] + p[3]*p[3];
		if(l > 0.07) return EF;
		else return p.getCoordSystem();
	}
	return -1;	//at this point apparently the point's coord system is invalid
}

/*
 * Metric in Eddington-Finkelstein coordinates
 */
 
SchwEFMetric::SchwEFMetric(int cS, SchwManifold* _m)
	: Metric(cS)
{ 
	m = _m;
}

SchwEFMetric::~SchwEFMetric()
{
}

double SchwEFMetric::_g(int i, int j, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	double M;
	M = m->getMass();
	
	double r = pos[coordR];
	double t = pos[coordTheta];
	
	if(i > j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k = i;
		i = j;
		j = k;
	}
	
	if(i == j)
	{
		if(i == coordU) return 1.0-2*M/r;
		if(i == coordR) return 0.0;
		if(i == coordTheta) return -r*r;
		if(i == coordPhi) return -r*r*sin(t)*sin(t);
	}
	if(i == coordU && j == coordR) return -1.0;
	
	return 0.0;
}

double SchwEFMetric::_invg(int i, int j, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	double M;
	M = m->getMass();
	
	double r = pos[coordR];
	double t = pos[coordTheta];
	
	if(i>j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k=i;
		i=j;
		j=k;
	}
	
	if(i == j)
	{
		if(i == coordR) return -1.0+2*M/r;
		if(i == coordTheta) return -1.0/r/r;
		if(i == coordPhi) return -1.0/(r*r*sin(t)*sin(t));
	}
	if(i == coordU && j == coordR) return -1;
	
	return 0.0;
}

double SchwEFMetric::_christoffel(int i, int j, int k, Point p)
{
	/*Point pos = m->convertPointTo(p, coordSystem);
	int n;
	
	double result = 0.0;
	
	for(n=0; n<4; n++)
		result += invg(i,n,p)*(dg(n,j,k,p) + dg(n,k,j,p) - dg(j,k,n,p));
	result *= 0.5;
	
	return result;*/
	
	Point pos = m->convertPointTo(p, coordSystem);
	double M;
	M = m->getMass();
	
	double r = pos[coordR];
	double t = pos[coordTheta];
	
	if(j>k)	//switch so that j<=k - it's symmetric anyway
	{
		int l;
		l=j;
		j=k;
		k=l;
	}
	
	switch(i)
	{
		case coordU:
			if(j == coordU && k == coordU) return M/r/r;
			if(j == coordTheta && k == coordTheta) return -r;
			if(j == coordPhi && k == coordPhi) return -r*sin(t)*sin(t);
			break;
		case coordR:
			if(j == coordU && k == coordU) return M*(1.0-2*M/r)/r/r;
			if(j == coordU && k == coordR) return -M/r/r;
			if(j == coordTheta && k == coordTheta) return -(r-2*M);
			if(j == coordPhi && k == coordPhi) return -(r-2*M)*sin(t)*sin(t);
			break;
		case coordTheta:
			if(j == coordR && k == coordTheta) return 1.0/r;
			if(j == coordPhi && k == coordPhi) return -sin(t)*cos(t);
			break;
		case coordPhi:
			if(j == coordR && k == coordPhi) return 1.0/r;
			if(j == coordTheta && k == coordPhi) return cos(t)/sin(t);
			break;
	}
	
	return 0.0;
}

/*
 * Metric in stereographic coordinates
 */
 
SchwNearPoleMetric::SchwNearPoleMetric(int cS, SchwManifold* _m)
	: Metric(cS)
{
	m = _m;
}

SchwNearPoleMetric::~SchwNearPoleMetric()
{
}


double SchwNearPoleMetric::_g(int i, int j, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	double M;
	M = m->getMass();
	
	if(i>j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k=i;
		i=j;
		j=k;
	}
	
	double r, x, y, alpha2;
	r = p[coordR];
	x = p[coordX];
	y = p[coordY];
	alpha2 = 4/(1.0+x*x+y*y)/(1.0+x*x+y*y);
	
	if(i == j)
	{
		if(i == coordU) return 1.0 - 2*M/r;
		if(i == coordX) return -alpha2*r*r;
		if(i == coordY) return -alpha2*r*r;
	}
	if(i == coordU && j == coordR) return -1.0;
	
	return 0.0;
}

double SchwNearPoleMetric::_invg(int i, int j, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	double M;
	M = m->getMass();
	
	if(i>j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k=i;
		i=j;
		j=k;
	}
	
	double r, x, y, alpha2;
	r = p[coordR];
	x = p[coordX];
	y = p[coordY];
	alpha2 = 4/(1.0+x*x+y*y)/(1.0+x*x+y*y);
	
	if(i == j)
	{
		if(i == coordR) return -(1.0 - 2*M/r);
		if(i == coordX || i == coordY) return -1.0/alpha2/r/r;
	}
	if(i == coordU && j == coordR) return -1.0;
	
	return 0.0;
}

double SchwNearPoleMetric::_christoffel(int i, int j, int k, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	int n;
	
	double result = 0.0;
	
	for(n=0; n<4; n++)
		result += invg(i,n,p)*(dg(n,j,k,p) + dg(n,k,j,p) - dg(j,k,n,p));
	result *= 0.5;
	
	return result;
}

