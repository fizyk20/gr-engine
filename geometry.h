#ifndef __GEOMETRY__
#define __GEOMETRY__

/*! \file geometry.h
 * \brief Basic geometry classes
 */

#ifndef NULL
#define NULL 0
#endif

class Manifold;

extern Manifold* gManifold;

/*! \class vector4
 * \brief A 4-dimensional tangent vector
 */
class vector4
{
	double v[4];
public:
	vector4();
	vector4(double, double, double, double);
	~vector4();
	
	vector4& operator=(const vector4&);
	vector4 operator+(const vector4);
	vector4 operator+=(const vector4);
	vector4 operator-(const vector4);
	vector4 operator-=(const vector4);
	vector4 operator*(const double);
	vector4 operator*=(const double);
	vector4 operator/(const double);
	vector4 operator/=(const double);
	friend vector4 operator*(const double, const vector4);
	
	double& operator[](int);
};

/*! \class Point
 * \brief A point on a manifold
 */
class Point
{
	double x[4];
	int coordSystem;
	static Manifold* m;
public:
	Point();
	Point(int cS, double x0=0.0, double x1=0.0, double x2=0.0, double x3=0.0);
	Point(int cS, double tx[]);
	Point(int cS, vector4 v);
	Point(const Point&);
	~Point();
	
	double& operator[](int);
	int getCoordSystem();
	vector4 toVector4();
	
	bool operator==(Point);
	bool operator!=(Point);
	
	static void setGlobalManifold(Manifold*);
};

/*! \class CoordinateConversion
 * \brief Basic class for conversions between coordinate systems
 */
class CoordinateConversion
{
public:
	CoordinateConversion();
	virtual ~CoordinateConversion() = 0;
	
	virtual Point convertPoint(Point) = 0;
	virtual double jacobian(int, int, Point) = 0;
	virtual double inv_jacobian(int, int, Point) = 0;
};

/*! \class IdentityConversion
 * \brief Class implementing the trivial conversion
 */
class IdentityConversion : public CoordinateConversion
{
public:
	IdentityConversion();
	~IdentityConversion();
	
	Point convertPoint(Point p);
	double jacobian(int i, int j, Point);
	double inv_jacobian(int i, int j, Point);
};

/*! \class Metric
 * \brief Class representing a metric on a manifold
 */
class Metric
{
	Point gCachePoints[4][4];
	Point invgCachePoints[4][4];
	Point gammaCachePoints[4][4][4];
	
	double gCache[4][4];
	double invgCache[4][4];
	double gammaCache[4][4][4];
protected:
	int coordSystem;
	double dg(int, int, int, Point);
	
	virtual double _g(int, int, Point) = 0;
	virtual double _invg(int, int, Point) = 0;
	virtual double _christoffel(int, int, int, Point) = 0;
public:
	Metric(int cS);
	virtual ~Metric() = 0;
	
	double g(int, int, Point);
	double invg(int, int, Point);
	double christoffel(int, int, int, Point);
	
	double g(vector4, vector4, Point);
	vector4 christoffel(vector4, vector4, Point);
};

/*! \class Manifold
 * \brief Class representing a manifold with different coordinate systems and metrics
 */
class Manifold
{
protected:
	int nCoordSystems;
	CoordinateConversion*** conversions;
	Metric** metrics;
public:
	Manifold();
	virtual ~Manifold() = 0;
	
	CoordinateConversion* getConversion(int, int);	//returns conversion from i to j
	Metric* getMetric(int);
	Point convertPointTo(Point, int);	//convert point to coordinates
	vector4 convertVectorTo(vector4, Point, int);	//convert vector to coordinates
	
	virtual int recommendCoordSystem(Point); //return the best coordinate system to use in a point
};

#endif

