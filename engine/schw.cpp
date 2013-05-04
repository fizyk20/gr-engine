#include "schw.h"
#include <math.h>

SchwManifold::SchwManifold(double _M)
{
	M = _M;
	
	gManifold = this;
	Point::setGlobalManifold(this);
	
	nCoordSystems = 3;
	
	metrics = new Metric*[nCoordSystems];
	metrics[EF] = new EFMetric(EF, this);
	metrics[NearPole0] = new NearPoleMetric(NearPole0, this);
	metrics[NearPolePi] = new NearPoleMetric(NearPolePi, this);
	
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
	case SchwManifold::EF:
		if(p[2] < 0.5) return SchwManifold::NearPole0;
		else if(p[2] > 2.642) return SchwManifold::NearPolePi;
		else return SchwManifold::EF;
	case SchwManifold::NearPole0:
	case SchwManifold::NearPolePi:
		l = p[2]*p[2] + p[3]*p[3];
		if(l > 0.07) return SchwManifold::EF;
		else return p.getCoordSystem();
	}
}

/*
 * Metric in Eddington-Finkelstein coordinates
 */
 
EFMetric::EFMetric(int cS, SchwManifold* _m)
	: Metric(cS)
{ 
	m = _m;
}

EFMetric::~EFMetric()
{
}

double EFMetric::_g(int i, int j, Point p)
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

double EFMetric::_invg(int i, int j, Point p)
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

double EFMetric::_christoffel(int i, int j, int k, Point p)
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
 
NearPoleMetric::NearPoleMetric(int cS, SchwManifold* _m)
	: Metric(cS)
{
	m = _m;
}

NearPoleMetric::~NearPoleMetric()
{
}


double NearPoleMetric::_g(int i, int j, Point p)
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

double NearPoleMetric::_invg(int i, int j, Point p)
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

double NearPoleMetric::_christoffel(int i, int j, int k, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	int n;
	
	double result = 0.0;
	
	for(n=0; n<4; n++)
		result += invg(i,n,p)*(dg(n,j,k,p) + dg(n,k,j,p) - dg(j,k,n,p));
	result *= 0.5;
	
	return result;
}

/*
 * Coordinate conversions
 */

EFToNearPole0::EFToNearPole0()
{
}

EFToNearPole0::~EFToNearPole0()
{
}

Point EFToNearPole0::convertPoint(Point p)
{
	if(p.getCoordSystem() != SchwManifold::EF) throw "EFToNearPole0: invalid coordinate system.";
	Point result(SchwManifold::NearPole0);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = tan(p[2]/2)*cos(p[3]);
	result[3] = tan(p[2]/2)*sin(p[3]);
	
	return result;
}

double EFToNearPole0::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::EF) throw "EFToNearPole0: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double t2 = p[2]/2;
		
	if(i == 2 && j == 2) return 2*cos(p[3])/(1.0+tan(t2)*tan(t2));
	if(i == 2 && j == 3) return 2*sin(p[3])/(1.0+tan(t2)*tan(t2));
	if(i == 3 && j == 2) return -sin(p[3])/tan(t2);
	if(i == 3 && j == 3) return cos(p[3])/tan(t2);
	
	return 0.0;
}

double EFToNearPole0::inv_jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::EF) throw "EFToNearPole0: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double t2 = p[2]/2;
	
	if(i == 2 && j == 2) return 0.5*cos(p[3])/cos(t2)/cos(t2);
	if(i == 2 && j == 3) return -tan(t2)*sin(p[3]);
	if(i == 3 && j == 2) return 0.5*sin(p[3])/cos(t2)/cos(t2);
	if(i == 3 && j == 3) return tan(t2)*cos(p[3]);
	
	return 0.0;
}

/*******/

EFToNearPolePi::EFToNearPolePi()
{
}

EFToNearPolePi::~EFToNearPolePi()
{
}

Point EFToNearPolePi::convertPoint(Point p)
{
	if(p.getCoordSystem() != SchwManifold::EF) throw "EFToNearPolePi: invalid coordinate system.";
	Point result(SchwManifold::NearPolePi);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = tan((M_PI-p[2])/2)*cos(p[3]);
	result[3] = tan((M_PI-p[2])/2)*sin(p[3]);
	
	return result;
}

double EFToNearPolePi::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::EF) throw "EFToNearPolePi: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double t2 = (M_PI-p[2])/2;
	
	if(i == 2 && j == 2) return -2*cos(p[3])/(1.0+tan(t2)*tan(t2));
	if(i == 2 && j == 3) return -2*sin(p[3])/(1.0+tan(t2)*tan(t2));
	if(i == 3 && j == 2) return -sin(p[3])/tan(t2);
	if(i == 3 && j == 3) return cos(p[3])/tan(t2);
	
	return 0.0;
}

double EFToNearPolePi::inv_jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::EF) throw "EFToNearPolePi: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double t2 = (M_PI-p[2])/2;
	
	if(i == 2 && j == 2) return -0.5*cos(p[3])/cos(t2)/cos(t2);
	if(i == 2 && j == 3) return -tan(t2)*sin(p[3]);
	if(i == 3 && j == 2) return -0.5*sin(p[3])/cos(t2)/cos(t2);
	if(i == 3 && j == 3) return tan(t2)*cos(p[3]);
	
	return 0.0;
}

/*******/

NearPole0ToEF::NearPole0ToEF()
{
}

NearPole0ToEF::~NearPole0ToEF()
{
}

Point NearPole0ToEF::convertPoint(Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	Point result(SchwManifold::EF);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = 2*atan(sqrt(p[2]*p[2]+p[3]*p[3]));
	result[3] = atan2(p[3],p[2]);
	
	return result;
}

double NearPole0ToEF::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
	
	double x,y;
	x = p[2];
	y = p[3];
	
	if(i == 2 && j == 2) return (1.0+x*x+y*y)/2*x/sqrt(x*x+y*y);
	if(i == 2 && j == 3) return -y;
	if(i == 3 && j == 2) return (1.0+x*x+y*y)/2*y/sqrt(x*x+y*y);
	if(i == 3 && j == 3) return x;
	
	return 0.0;
}

double NearPole0ToEF::inv_jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double x,y;
	x = p[2];
	y = p[3];
	
	if(i == 2 && j == 2) return 2.0/(1.0+x*x+y*y)*x/sqrt(x*x+y*y);
	if(i == 2 && j == 3) return 2.0/(1.0+x*x+y*y)*y/sqrt(x*x+y*y);
	if(i == 3 && j == 2) return -y/(x*x+y*y);
	if(i == 3 && j == 3) return x/(x*x+y*y);
	
	return 0.0;
}

/*******/

NearPolePiToEF::NearPolePiToEF()
{
}

NearPolePiToEF::~NearPolePiToEF()
{
}

Point NearPolePiToEF::convertPoint(Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	Point result(SchwManifold::EF);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = M_PI - 2*atan(sqrt(p[2]*p[2]+p[3]*p[3]));
	result[3] = atan2(p[3],p[2]);
	
	return result;
}

double NearPolePiToEF::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double x,y;
	x = p[2];
	y = p[3];
	
	if(i == 2 && j == 2) return -(1.0+x*x+y*y)/2*x/sqrt(x*x+y*y);
	if(i == 2 && j == 3) return -y;
	if(i == 3 && j == 2) return -(1.0+x*x+y*y)/2*y/sqrt(x*x+y*y);
	if(i == 3 && j == 3) return x;
	
	return 0.0;
}

double NearPolePiToEF::inv_jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double x,y;
	x = p[2];
	y = p[3];
	
	if(i == 2 && j == 2) return -2.0/(1.0+x*x+y*y)*x/sqrt(x*x+y*y);
	if(i == 2 && j == 3) return -2.0/(1.0+x*x+y*y)*y/sqrt(x*x+y*y);
	if(i == 3 && j == 2) return -y/(x*x+y*y);
	if(i == 3 && j == 3) return x/(x*x+y*y);
	
	return 0.0;
}

/*******/

NearPoleToNearPole::NearPoleToNearPole()
{
}

NearPoleToNearPole::~NearPoleToNearPole()
{
}

Point NearPoleToNearPole::convertPoint(Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPole0 && p.getCoordSystem() != SchwManifold::NearPolePi) 
		throw "NearPoleToNearPole: invalid coordinate system.";
	
	Point result(p.getCoordSystem() == SchwManifold::NearPole0 ? SchwManifold::NearPolePi : SchwManifold::NearPole0);
	
	double l = p[2]*p[2] + p[3]*p[3];
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = p[2]/l;
	result[3] = p[3]/l;
	
	return result;
}

double NearPoleToNearPole::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPole0 && p.getCoordSystem() != SchwManifold::NearPolePi) 
		throw "NearPoleToNearPole: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double x,y;
	x = p[2];
	y = p[3];
	double l = x*x+y*y;
	
	if(i == 2 && j == 2) return (y*y-x*x)/l;
	if(i == 2 && j == 3) return -2*x*y/l;
	if(i == 3 && j == 2) return -2*x*y/l;
	if(i == 3 && j == 3) return (x*x-y*y)/l;
	
	return 0.0;
}

double NearPoleToNearPole::inv_jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != SchwManifold::NearPole0 && p.getCoordSystem() != SchwManifold::NearPolePi) 
		throw "NearPoleToNearPole: invalid coordinate system.";
	
	if(i == 0 || i == 1)
	{
		if(i == j) return 1.0;
		else return 0.0;
	}
		
	double x,y;
	x = p[2];
	y = p[3];
	double l = x*x+y*y;
	
	if(i == 2 && j == 2) return (y*y-x*x)/l;
	if(i == 2 && j == 3) return -2*x*y/l;
	if(i == 3 && j == 2) return -2*x*y/l;
	if(i == 3 && j == 3) return (x*x-y*y)/l;
	
	return 0.0;
}

