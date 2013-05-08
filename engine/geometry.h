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
 *
 * The class represents a 4-vector. It supports basic operation like addition, subtraction and multiplication and division by a scalar.
 */
class vector4
{
	double v[4];	///< The components of the vector.
public:
	//! Default onstructor
	/*! Constructs a vector (0, 0, 0, 0)
	 */
	vector4();
	//! Constructor
	/*! \param v1 First component
	 *  \param v2 Second component
	 *  \param v3 Third component
	 *  \param v4 Fourth component
	 */
	vector4(double v1, double v2, double v3, double v4);
	//! Destructor
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
	
	//! Index operator
	/*! \param i Number of the component to retrieve
	    \return Reference to the component, which allows for modification
	 */
	double& operator[](int i);
};

/*! \class Point
 * \brief A point on a manifold
 *
 * The class defines a point on the manifold using 4 coordinates and the number of the coordinate system being used. This allows for two points from different coordinate systems to be compared.
 */
class Point
{
	double x[4];
	int coordSystem;
	static Manifold* m; ///< The manifold used by the program
public:
	//! Default constructor - initializes the point with an invalid coordinate system
	Point();
	//! Constructor
	/*! \param cS Coordinate system being used
	 *  \param x0 First coordinate
	 *  \param x1 Second coordinate
	 *  \param x2 Third coordinate
	 *  \param x3 Fourth coordinate
	 */
	Point(int cS, double x0=0.0, double x1=0.0, double x2=0.0, double x3=0.0);
	//! Constructor
	/*! \param cS Coordinate system being used
	 *  \param tx Array of coordinates
	 */
	Point(int cS, double tx[]);
	//! Constructor
	/*! \param cS Coordinate system being used
	 *  \param v A vector4 containing the coordinates
	 */
	Point(int cS, vector4 v);
	//! Copy constructor
	Point(const Point&);
	//! Destructor
	~Point();
	
	//! Index operator
	/*! \param i Number of the coordinate to retrieve
	    \return Reference to the coordinate, which allows for modification
	 */
	double& operator[](int i);
	//! Function returning the coordinate system
	/* \return The number of the coordinate system in use
	 */
	int getCoordSystem();
	//! Conversion to vector4
	/* \return A vector4 with components equal to the coordinates of the point
	 */
	vector4 toVector4();
	
	bool operator==(Point);
	bool operator!=(Point);
	
	//! Static method setting the global manifold variable
	/* This method sets the manifold being used in the whole program
	 * \param _m The pointer to the manifold object
	 */
	static void setGlobalManifold(Manifold* _m);
};

/*! \class CoordinateConversion
 * \brief Basic class for conversions between coordinate systems
 *
 * This class provides interface for converting points and vectors between different coordinate systems defined on a manifold.
 */
class CoordinateConversion
{
public:
	//! Constructor
	CoordinateConversion();
	//! Destructor
	virtual ~CoordinateConversion();
	
	//! Conversion of a point
	/*! This method converts a single point to a different coordinate system.
	 *  \param p Point to be converted
	 *  \return Converted point
	 */
	virtual Point convertPoint(Point p) = 0;
	//! Jacobian of the conversion
	/*! This method returns the jacobian (matrix of partial derivatives) of the transformation between coordinate systems.
	 *  \param i The number of the row
	 *  \param j The number of the column
	 *  \param p The point at which the jacobian is calculated
	 *  \return dx^i/dy^j(p), where x^i - old coordinates, y^j - new coordinates
	 */
	virtual double jacobian(int i, int j, Point p) = 0;
	//! Inverse jacobian of the conversion
	/*! This method returns the inverse jacobian (matrix of partial derivatives) of the transformation between coordinate systems.
	 *  \param i The number of the row
	 *  \param j The number of the column
	 *  \param p The point at which the jacobian is calculated
	 *  \return dy^i/dx^j(p), where y^i - new coordinates, x^j - old coordinates
	 */
	virtual double inv_jacobian(int i, int j, Point p) = 0;
};

/*! \class IdentityConversion
 * \brief Class implementing the trivial conversion
 *
 * This is the predefined class implementing the conversion given by: new_coord_i = old_coord_i
 */
class IdentityConversion : public CoordinateConversion
{
public:
	//! Constructor
	IdentityConversion();
	//! Destructor
	~IdentityConversion();
	
	//! Conversion of a point
	/*! This method converts a single point to a different coordinate system.
	 *  \param p Point to be converted
	 *  \return Converted point (the same as p)
	 */
	Point convertPoint(Point p);
	//! Jacobian of the conversion
	/*! This method returns the jacobian (matrix of partial derivatives) of the transformation between coordinate systems (equal to the identity matrix)
	 *  \param i The number of the row
	 *  \param j The number of the column
	 *  \param p The point at which the jacobian is calculated
	 *  \return delta_ij (Kronecker delta)
	 */
	double jacobian(int i, int j, Point);
	//! Inverse jacobian of the conversion
	/*! This method returns the inverse jacobian (matrix of partial derivatives) of the transformation between coordinate systems (equal to the identity matrix)
	 *  \param i The number of the row
	 *  \param j The number of the column
	 *  \param p The point at which the jacobian is calculated
	 *  \return delta_ij (Kronecker delta)
	 */
	double inv_jacobian(int i, int j, Point);
};

/*! \class Metric
 * \brief Class representing a metric on a manifold
 *
 * This class represents the metric tensor on a manifold, expressed in some coordinate system.
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
	//! Partial derivative of the metric
	/*! Returns partial derivative of a component of the metric with respect to some coordinate.
	 *  \param i First index of the component
	 *  \param j Second index of the component
	 *  \param k The index of the coordinate, with respect to which the derivative is calculated
	 *  \param p The point at which the derivative is evaluated
	 *  \return dg_ij/dx^k(p)
	 */
	double dg(int i, int j, int k, Point p);
	
	//! Component of the metric
	/*! Function calculating a component of the metric.
	 *  \param i First index of the component
	 *  \param j Second index of the component
	 *  \param p The point at which the component is evaluated
	 *  \return g_ij(p)
	 */
	virtual double _g(int i, int j, Point p) = 0;
	//! Component of the inverse metric
	/*! Function calculating a component of the inverse metric.
	 *  \param i First index of the component
	 *  \param j Second index of the component
	 *  \param p The point at which the component is evaluated
	 *  \return g^ij(p)
	 */
	virtual double _invg(int i, int j, Point p) = 0;
	//! Component of the Christoffel symbol
	/*! Function calculating a component of the Christoffel symbol.
	 *  \param i First index of the component
	 *  \param j Second index of the component
	 *  \param k Third index of the component
	 *  \param p The point at which the component is evaluated
	 *  \return Gamma^i_jk(p)
	 */
	virtual double _christoffel(int i, int j, int k, Point p) = 0;
public:
    //! Constructor
    /*! \param cS Coordinate system.
     */
	Metric(int cS);
	//! Destructor
	virtual ~Metric();
	
	//! Component of the metric
	/*! Function returning a component of the metric
	 *  Uses internal caching.
	 *  \param i First index of the component
	 *  \param j Second index of the component
	 *  \param p The point at which the component is evaluated
	 *  \return g_ij(p)
	 */
	double g(int i, int j, Point p);
	//! Component of the inverse metric
	/*! Function returning a component of the inverse metric
	 *  Uses internal caching.
	 *  \param i First index of the component
	 *  \param j Second index of the component
	 *  \param p The point at which the component is evaluated
	 *  \return g^ij(p)
	 */
	double invg(int i, int j, Point p);
	//! Component of the Christoffel symbol
	/*! Function returning a component of the Christoffel symbol
	 *  Uses internal caching.
	 *  \param i First index of the component
	 *  \param j Second index of the component
	 *  \param k Third index of the component
	 *  \param p The point at which the component is evaluated
	 *  \return Gamma^i_jk(p)
	 */
	double christoffel(int i, int j, int k, Point p);
	
	//! Dot product
	/*! Function returning the dot product of two vectors.
	 *  \param u First vector
	 *  \param v Second vector
	 *  \param p The point at which the dot product is evaluated
	 *  \return g_ij(p) u^i v^j
	 */
	double g(vector4 u, vector4 v, Point p);
	
	//! Christoffel symbol acting on two vectors
	/*! Function returning result of contraction of the Christoffel symbol with two vectors.
	 *  \param u First vector
	 *  \param v Second vector
	 *  \param p The point at which the dot product is evaluated
	 *  \return Gamma^i_jk(p) u^j v^k
	 */
	vector4 christoffel(vector4 u, vector4 v, Point p);
};

/*! \class Manifold
 * \brief Class representing a manifold with different coordinate systems and metrics
 *
 *  This class represents a manifold with multiple coordinate systems, conversions between them and a metric, which can be expressed in either of the systems.
 */
class Manifold
{
protected:
	int nCoordSystems;
	CoordinateConversion*** conversions;
	Metric** metrics;
public:
	//! Constructor
	Manifold();
	//! Destructor
	virtual ~Manifold();
	
	//! Get a coordinate conversion
	/*! \param i Number of the source coordinate system
	 *  \param j Number of the target coordinate system
	 *  \return Pointer to the conversion object
	 */
	CoordinateConversion* getConversion(int i, int j);	//returns conversion from i to j
	//! Get a metric
	/*! \param i Number of the coordinate system
	 *  \return Pointer to the metric expressed in the chosen system
	 */
	Metric* getMetric(int i);
	//! Convert a point to another coordinate system
	/*! \param p The point to be converted
	 *  \param system The target coordinate system
	 *  \return Converted point
	 */
	Point convertPointTo(Point p, int system);	//convert point to coordinates
	//! Convert a vector to another coordinate system
	/*! \param v The vector to be converted
	 *  \param p The point at which the vector is defined
	 *  \param system The target coordinate system
	 *  \return Converted vector
	 */
	vector4 convertVectorTo(vector4 v, Point p, int system);	//convert vector to coordinates
	
	//! Get the best coordinate system to use at a point
	/*! \param p The point
	 *  \return The recommended coordinate system to be used
	 */
	virtual int recommendCoordSystem(Point p);
};

#endif

