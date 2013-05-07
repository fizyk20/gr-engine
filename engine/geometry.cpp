#include "geometry.h"

/*
Point
*/

Manifold* Point::m;

Point::Point()
{
	coordSystem = -1;
}

Point::Point(int cS, double x0, double x1, double x2, double x3)
{
	coordSystem = cS;
	x[0] = x0;
	x[1] = x1;
	x[2] = x2;
	x[3] = x3;
}

Point::Point(int cS, double tx[])
{
	coordSystem = cS;
	int i;
	for(i=0; i<4; i++)
		x[i] = tx[i];
}

Point::Point(int cS, vector4 v)
{
	coordSystem = cS;
	int i;
	for(i=0; i<4; i++)
		x[i] = v[i];
}

Point::Point(const Point& p)
{
	coordSystem = p.coordSystem;
	int i;
	for(i=0; i<4; i++)
		x[i] = p.x[i];
}

Point::~Point()
{
}

double& Point::operator[](int i)
{
	if(i<0 || i>3) throw "Point: Index out of bounds.";
	return x[i];
}

int Point::getCoordSystem()
{
	return coordSystem;
}

vector4 Point::toVector4()
{
	return vector4(x[0], x[1], x[2], x[3]);
}

bool Point::operator==(Point arg)
{
	if(coordSystem == -1) return false;
	if(arg.coordSystem == -1) return false;
	
	if(arg.coordSystem != coordSystem)
		arg = m->getConversion(arg.coordSystem, coordSystem)->convertPoint(arg);
	
	for(int i=0; i<4; i++)
		if(x[i] != arg.x[i]) return false;
	return true;
}

bool Point::operator!=(Point arg)
{
	return !((*this) == arg);
}

void Point::setGlobalManifold(Manifold* _m)
{
	m = _m;
}

/*
vector4
*/

vector4::vector4()
{
	v[0] = v[1] = v[2] = v[3] = 0.0;
}

vector4::vector4(double v1, double v2, double v3, double v4)
{
	v[0] = v1;
	v[1] = v2;
	v[2] = v3;
	v[3] = v4;
}

vector4::~vector4()
{
}

vector4& vector4::operator=(const vector4& arg)
{
	v[0] = arg.v[0];
	v[1] = arg.v[1];
	v[2] = arg.v[2];
	v[3] = arg.v[3];
	
	return (*this);
}

vector4 vector4::operator+(const vector4 arg)
{
	vector4 result(v[0]+arg.v[0], v[1]+arg.v[1], v[2]+arg.v[2], v[3]+arg.v[3]);
	return result;
}

vector4 vector4::operator+=(const vector4 arg)
{
	v[0] += arg.v[0];
	v[1] += arg.v[1];
	v[2] += arg.v[2];
	v[3] += arg.v[3];
	
	return (*this);
}

vector4 vector4::operator-(const vector4 arg)
{
	vector4 result(v[0]-arg.v[0], v[1]-arg.v[1], v[2]-arg.v[2], v[3]-arg.v[3]);
	return result;
}

vector4 vector4::operator-=(const vector4 arg)
{
	v[0] -= arg.v[0];
	v[1] -= arg.v[1];
	v[2] -= arg.v[2];
	v[3] -= arg.v[3];
	
	return (*this);
}

vector4 vector4::operator*(const double arg)
{
	vector4 result(v[0]*arg, v[1]*arg, v[2]*arg, v[3]*arg);
	return result;
}

vector4 vector4::operator*=(const double arg)
{
	v[0] *= arg;
	v[1] *= arg;
	v[2] *= arg;
	v[3] *= arg;
	
	return (*this);
}

vector4 vector4::operator/(const double arg)
{
	vector4 result(v[0]/arg, v[1]/arg, v[2]/arg, v[3]/arg);
	return result;
}

vector4 vector4::operator/=(const double arg)
{
	v[0] /= arg;
	v[1] /= arg;
	v[2] /= arg;
	v[3] /= arg;
	
	return (*this);
}

vector4 operator*(const double a, const vector4 v)
{
	vector4 res(v.v[0]*a, v.v[1]*a, v.v[2]*a, v.v[3]*a);
	return res;
}

double& vector4::operator[](int i)
{
	if(i<0) return v[0];
	if(i>3) return v[3];
	return v[i];
}

/*
 * CoordinateConversion
 */
 
CoordinateConversion::CoordinateConversion()
{
}

CoordinateConversion::~CoordinateConversion()
{
}

/*
 * IdentityConversion
 */
 
IdentityConversion::IdentityConversion()
{
}
 
IdentityConversion::~IdentityConversion()
{
}

Point IdentityConversion::convertPoint(Point p)
{ 
	return p;
}

double IdentityConversion::jacobian(int i, int j, Point)
{
	return (i == j) ? 1.0 : 0.0;
}

double IdentityConversion::inv_jacobian(int i, int j, Point)
{ 
	return (i == j) ? 1.0 : 0.0;
}

/*
 * Metric
 */
 
Metric::Metric(int cS)
{
	coordSystem = cS;
}

Metric::~Metric()
{
}

double Metric::dg(int i, int j, int k, Point p)
{
	double h = 0.0001;
	double df = 0.0;
	
	p[k] += 2*h;
	df -= g(i,j,p);
	p[k] -= h;
	df += 8*g(i,j,p);
	p[k] -= 2*h;
	df -= 8*g(i,j,p);
	p[k] -= h;
	df += g(i,j,p);
	
	return df/(12*h);
}

double Metric::g(vector4 u, vector4 v, Point p)
{
	int i,j;
	double sum = 0.0;
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			sum += u[i]*v[j]*g(i,j,p);
	return sum;
}

vector4 Metric::christoffel(vector4 u, vector4 v, Point p)
{
	int i, j, k;
	vector4 result;
	
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			for(k=0; k<4; k++)
				result[i] += u[j]*v[k]*christoffel(i,j,k,p);
	return result;
}

double Metric::g(int i, int j, Point p)
{
	if(i > j)
	{
		int k;
		k = i;
		i = j;
		j = k;
	}
	
	if(gCachePoints[i][j] != p)
	{
		gCachePoints[i][j] = p;
		gCache[i][j] = _g(i, j, p);
	}
	return gCache[i][j];
}

double Metric::invg(int i, int j, Point p)
{
	if(i > j)
	{
		int k;
		k = i;
		i = j;
		j = k;
	}
	
	if(invgCachePoints[i][j] != p)
	{
		invgCachePoints[i][j] = p;
		invgCache[i][j] = _invg(i, j, p);
	}
	return invgCache[i][j];
}

double Metric::christoffel(int i, int j, int k, Point p)
{
	if(j > k)
	{
		int l;
		l = j;
		j = k;
		k = l;
	}
	
	if(gammaCachePoints[i][j][k] != p)
	{
		gammaCachePoints[i][j][k] = p;
		gammaCache[i][j][k] = _christoffel(i, j, k, p);
	}
	return gammaCache[i][j][k];
}

/*
Manifold
*/

Manifold* gManifold = NULL;

Manifold::Manifold()
{
}

Manifold::~Manifold()
{
}

CoordinateConversion* Manifold::getConversion(int i, int j)
{
	if(i < 0 || i >= nCoordSystems || j < 0 || j >= nCoordSystems) throw "Manifold: Index out of bounds.";
	return conversions[i][j];
}

Metric* Manifold::getMetric(int i)
{
	if(i < 0 || i >= nCoordSystems) throw "Manifold: Index out of bounds.";
	return metrics[i];
}

Point Manifold::convertPointTo(Point p, int system)
{
	if(system < 0 || system >= nCoordSystems) throw "Manifold: Index out of bounds.";
	Point result(system);
	
	result = conversions[p.getCoordSystem()][system]->convertPoint(p);
	
	return result;
}

vector4 Manifold::convertVectorTo(vector4 v, Point p, int system)
{
	if(system < 0 || system >= nCoordSystems) throw "Manifold: Index out of bounds.";
	vector4 result;
	
	CoordinateConversion* c = conversions[p.getCoordSystem()][system];
	int i,j;
	
	for(i=0; i<4; i++)
	{
		result[i] = 0;
		for(j=0; j<4; j++)
			result[i] += c->inv_jacobian(i,j,p)*v[j];
	}
	
	return result;
}

int Manifold::recommendCoordSystem(Point p)
{
	return p.getCoordSystem(); 	//trivial implementation
}
