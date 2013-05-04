#include "kerr.h"
#include <math.h>

KerrManifold::KerrManifold(double _M, double _a)
{
	M = _M;
	a = _a;
	
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

KerrManifold::~KerrManifold()
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
	
double KerrManifold::getMass()
{
	return M;
}

double KerrManifold::getAngMomentum()
{
	return a;
}

void KerrManifold::setMass(double _M)
{
	M = _M;
}

void KerrManifold::setAngMomentum(double _a)
{
	a = _a;
}
	
int KerrManifold::recommendCoordSystem(Point p)
{
	double l;
	switch(p.getCoordSystem())
	{
	case KerrManifold::EF:
		if(p[2] < 0.5) return KerrManifold::NearPole0;
		else if(p[2] > 2.642) return KerrManifold::NearPolePi;
		else return KerrManifold::EF;
	case KerrManifold::NearPole0:
	case KerrManifold::NearPolePi:
		l = p[2]*p[2] + p[3]*p[3];
		if(l > 0.07) return KerrManifold::EF;
		else return p.getCoordSystem();
	default:
		return 0;
	}
}

/*
 * Metric in Eddington-Finkelstein coordinates
 */
 
EFMetric::EFMetric(int cS, KerrManifold* _m)
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
	double M,a;
	M = m->getMass();
	a = m->getAngMomentum();
	
	double r = pos[coordR];
	double t = pos[coordTheta];
	
	if(i > j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k = i;
		i = j;
		j = k;
	}
	
	double rho2 = r*r + a*a*cos(t)*cos(t);
	
	if(i == j)
	{
		if(i == coordU) return 1.0-2*M*pos[1]/rho2;
		if(i == coordR) return 0.0;
		if(i == coordTheta) return -rho2;
		if(i == coordPhi) return -(r*r+a*a+2*M*r*a*a*sin(t)*sin(t)/rho2)*sin(t)*sin(t);
	}
	if(i == coordU && j == coordR) return -1.0;
	if(i == coordU && j == coordPhi) return 2*M*r*a*sin(t)*sin(t)/rho2;
	if(i == coordR && j == coordPhi) return a*sin(t)*sin(t);
	
	return 0.0;
}

double EFMetric::_invg(int i, int j, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	double M,a;
	M = m->getMass();
	a = m->getAngMomentum();
	
	double r = pos[coordR];
	double t = pos[coordTheta];
	
	double rho2 = r*r + a*a*cos(t)*cos(t);
	
	if(i>j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k=i;
		i=j;
		j=k;
	}
	
	if(i == j)
	{
		if(i == coordU) return -a*a/rho2*sin(t)*sin(t);
		if(i == coordR)
		{
			double delta = r*r-2*M*r+a*a;
			return -delta/rho2;
		}
		if(i == coordTheta) return -1.0/rho2;
		if(i == coordPhi) return -1.0/(rho2*sin(t)*sin(t));
	}
	if(i == coordU && j == coordR) return -(r*r+a*a)/rho2;
	if((i == coordU && j == coordPhi) || (i == coordR && j == coordPhi)) return -a/rho2;
	
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
	double M,a;
	M = m->getMass();
	a = m->getAngMomentum();
	
	double r = pos[coordR];
	double t = pos[coordTheta];
	
	double rho2 = r*r + a*a*cos(t)*cos(t);
	double rho2x = r*r - a*a*cos(t)*cos(t);
	double rho4 = rho2*rho2;
	double rho6 = rho2*rho4;
	double delta = r*r-2*M*r+a*a;
	
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
			if(j == coordU && k == coordU) return M*(r*r+a*a)*rho2x/rho6;
			if(j == coordU && k == coordTheta) return -M*r*a*a*sin(2*t)/rho4;
			if(j == coordU && k == coordPhi) return -M*a*(r*r+a*a)*rho2x*sin(t)*sin(t)/rho6;
			if(j == coordR && k == coordTheta) return -a*a*sin(t)*cos(t)/rho2;
			if(j == coordR && k == coordPhi) return a*r*sin(t)*sin(t)/rho2;
			if(j == coordTheta && k == coordTheta) return -(r*r+a*a)/rho2*r;
			if(j == coordTheta && k == coordPhi) return 2*M*r*a*a*a*sin(t)*sin(t)*sin(t)*cos(t)/rho4;
			if(j == coordPhi && k == coordPhi) return (r*r+a*a)/rho2*(M*rho2x*a*a*sin(t)*sin(t)*sin(t)*sin(t)/rho4-r*sin(t)*sin(t));
			break;
		case coordR:
			if(j == coordU && k == coordU) return M*rho2x*delta/rho6;
			if(j == coordU && k == coordR) return -M*rho2x/rho4;
			if(j == coordU && k == coordPhi) return -M*rho2x*delta*a*sin(t)*sin(t)/rho6;
			if(j == coordR && k == coordTheta) return -a*a*sin(t)*cos(t)/rho2;
			if(j == coordR && k == coordPhi) return (r*rho2+M*rho2x)*a*sin(t)*sin(t)/rho4;
			if(j == coordTheta && k == coordTheta) return -delta/rho2*r;
			if(j == coordPhi && k == coordPhi) return delta*sin(t)*sin(t)*(M*a*a*rho2x*sin(t)*sin(t)-r*rho4)/rho6;
			break;
		case coordTheta:
			if(j == coordU && k == coordU) return -2*M*r*a*a*sin(t)*cos(t)/rho6;
			if(j == coordU && k == coordPhi) return 2*M*r*a*(r*r+a*a)*sin(t)*cos(t)/rho6;
			if(j == coordR && k == coordTheta) return r/rho2;
			if(j == coordR && k == coordPhi) return a*sin(t)*cos(t)/rho2;
			if(j == coordTheta && k == coordTheta) return -a*a*sin(t)*cos(t)/rho2;
			if(j == coordPhi && k == coordPhi) return -sin(t)*cos(t)*(rho4*(r*r+a*a)+2*M*r*a*a*sin(t)*sin(t)*(r*r+a*a+rho2))/rho6;
			break;
		case coordPhi:
			if(j == coordU && k == coordU) return M*a*rho2x/rho6;
			if(j == coordU && k == coordTheta) return -2*M*r*a*cos(t)/(rho4*sin(t));
			if(j == coordU && k == coordPhi) return -M*a*a*rho2x*sin(t)*sin(t)/rho6;
			if(j == coordR && k == coordTheta) return -a/rho2*cos(t)/sin(t);
			if(j == coordR && k == coordPhi) return r/rho2;
			if(j == coordTheta && k == coordTheta) return -a*r/rho2;
			if(j == coordTheta && k == coordPhi) return cos(t)/sin(t)*(1+2*M*r*a*a*sin(t)*sin(t)/rho4);
			if(j == coordPhi && k == coordPhi) return a*sin(t)*sin(t)*(M*a*a*rho2x*sin(t)*sin(t)-r*rho4)/rho6;
			break;
	}
	
	return 0.0;
}

/*
 * Metric in stereographic coordinates
 */
 
NearPoleMetric::NearPoleMetric(int cS, KerrManifold* _m)
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
	double M,a;
	M = m->getMass();
	a = m->getAngMomentum();
	
	if(i>j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k=i;
		i=j;
		j=k;
	}
	
	double rho2, r, x, y, alpha2;
	r = p[coordR];
	x = p[coordX];
	y = p[coordY];
	rho2 = r*r + a*a*(1.0-x*x-y*y)*(1.0-x*x-y*y)/(1.0+x*x+y*y)/(1.0+x*x+y*y);
	alpha2 = 4/(1.0+x*x+y*y)/(1.0+x*x+y*y);
	
	if(i == j)
	{
		if(i == coordU) return 1.0 - 2*M*r/rho2;
		if(i == coordX) return -alpha2*(r*r + a*a - alpha2*a*a*(x*x - 2*M*r*y*y/rho2));
		if(i == coordY) return -alpha2*(r*r + a*a - alpha2*a*a*(y*y - 2*M*r*x*x/rho2));
	}
	if(i == coordU && j == coordR) return -1.0;
	if(i == coordU && j == coordX) return -2*M*r*a*y*alpha2/rho2;
	if(i == coordU && j == coordY) return 2*M*r*a*x*alpha2/rho2;
	if(i == coordR && j == coordX) return -alpha2*a*y;
	if(i == coordR && j == coordY) return alpha2*a*x;
	if(i == coordX && j == coordY) return x*y*a*a*alpha2*alpha2*(1.0+2*M*r/rho2);
	
	return 0.0;
}

double NearPoleMetric::_invg(int i, int j, Point p)
{
	Point pos = m->convertPointTo(p, coordSystem);
	double M,a;
	M = m->getMass();
	a = m->getAngMomentum();
	
	if(i>j)	//switch so that i<=j - it's symmetric anyway
	{
		int k;
		k=i;
		i=j;
		j=k;
	}
	
	double rho2, r, x, y, alpha2;
	r = p[coordR];
	x = p[coordX];
	y = p[coordY];
	rho2 = r*r + a*a*(1.0-x*x-y*y)*(1.0-x*x-y*y)/(1.0+x*x+y*y)/(1.0+x*x+y*y);
	alpha2 = 4/(1.0+x*x+y*y)/(1.0+x*x+y*y);
	
	if(i == j)
	{
		if(i == coordU) return -alpha2*a*a*(x*x+y*y)/rho2;
		if(i == coordR) return -(r*r + a*a - 2*M*r)/rho2;
		if(i == coordX || i == coordY) return -1.0/alpha2/rho2;
	}
	if(i == coordU && j == coordR) return -(r*r+a*a)/rho2;
	if(i == coordU && j == coordX) return a*y/rho2;
	if(i == coordU && j == coordY) return -a*x/rho2;
	if(i == coordR && j == coordX) return a*y/rho2;
	if(i == coordR && j == coordY) return -a*x/rho2;
	
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
	if(p.getCoordSystem() != KerrManifold::EF) throw "EFToNearPole0: invalid coordinate system.";
	Point result(KerrManifold::NearPole0);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = tan(p[2]/2)*cos(p[3]);
	result[3] = tan(p[2]/2)*sin(p[3]);
	
	return result;
}

double EFToNearPole0::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != KerrManifold::EF) throw "EFToNearPole0: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::EF) throw "EFToNearPole0: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::EF) throw "EFToNearPolePi: invalid coordinate system.";
	Point result(KerrManifold::NearPolePi);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = tan((M_PI-p[2])/2)*cos(p[3]);
	result[3] = tan((M_PI-p[2])/2)*sin(p[3]);
	
	return result;
}

double EFToNearPolePi::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != KerrManifold::EF) throw "EFToNearPolePi: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::EF) throw "EFToNearPolePi: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	Point result(KerrManifold::EF);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = 2*atan(sqrt(p[2]*p[2]+p[3]*p[3]));
	result[3] = atan2(p[3],p[2]);
	
	return result;
}

double NearPole0ToEF::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != KerrManifold::NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	Point result(KerrManifold::EF);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = M_PI - 2*atan(sqrt(p[2]*p[2]+p[3]*p[3]));
	result[3] = atan2(p[3],p[2]);
	
	return result;
}

double NearPolePiToEF::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != KerrManifold::NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::NearPole0 && p.getCoordSystem() != KerrManifold::NearPolePi) 
		throw "NearPole0ToNearPolePi: invalid coordinate system.";
	
	Point result(p.getCoordSystem() == KerrManifold::NearPole0 ? KerrManifold::NearPolePi : KerrManifold::NearPole0);
	
	double l = p[2]*p[2] + p[3]*p[3];
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = p[2]/l;
	result[3] = p[3]/l;
	
	return result;
}

double NearPoleToNearPole::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != KerrManifold::NearPole0 && p.getCoordSystem() != KerrManifold::NearPolePi) 
		throw "NearPole0ToNearPolePi: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != KerrManifold::NearPole0 && p.getCoordSystem() != KerrManifold::NearPolePi) 
		throw "NearPole0ToNearPolePi: invalid coordinate system.";
	
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

