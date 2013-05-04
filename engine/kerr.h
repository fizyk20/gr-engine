#ifndef __KERR_H__
#define __KERR_H__

/*! \file kerr.h
 * \brief Implementation of the Kerr metric
 */

#include "geometry.h"

/*! \class KerrManifold
 * \brief Class representing a Kerr manifold
 */
class KerrManifold : public Manifold
{
	double M, a;
public:
	enum { EF = 0, NearPole0 = 1, NearPolePi = 2 };
	
	KerrManifold(double, double);
	~KerrManifold();
	
	double getMass();
	double getAngMomentum();
	void setMass(double _M);
	void setAngMomentum(double _a);
	
	int recommendCoordSystem(Point);
};

/*! \class EFMetric
 * \brief The Schwarzschild metric in Eddington-Finkelstein coordinates
 */
class EFMetric : public Metric
{
	KerrManifold* m;
protected:
	double _g(int, int, Point);
	double _invg(int, int, Point);
	double _christoffel(int, int, int, Point);
	
public:
	enum { coordU = 0, coordR = 1, coordTheta = 2, coordPhi = 3 };
	
	EFMetric(int cS, KerrManifold* _m);
	~EFMetric();
};

/*! \class NearPoleMetric
 * \brief The Schwarzschild metric valid near spherical poles - in "stereographic" coordinates
 */
class NearPoleMetric : public Metric
{
	KerrManifold* m;
protected:
	double _g(int, int, Point);
	double _invg(int, int, Point);
	double _christoffel(int, int, int, Point);
	
public:
	enum { coordU = 0, coordR = 1, coordX = 2, coordY = 3 };
	
	NearPoleMetric(int cS, KerrManifold* _m);
	~NearPoleMetric();
};

/*! \class EFToNearPole0
 * \brief Coordinate conversion from EF coordinates to stereographic near pole theta=0
 */
class EFToNearPole0 : public CoordinateConversion
{
public:
	EFToNearPole0();
	~EFToNearPole0();
	
	Point convertPoint(Point);
	double jacobian(int i, int j, Point);	//d(old_i)/d(new_j)
	double inv_jacobian(int i, int j, Point);	//d(new_i)/d(old_j)
};

/*! \class NearPole0ToEF
 * \brief Coordinate conversion from stereographic coordinates near pole theta=0 to EF
 */
class NearPole0ToEF : public CoordinateConversion
{
public:
	NearPole0ToEF();
	~NearPole0ToEF();
	
	Point convertPoint(Point);
	double jacobian(int, int, Point);
	double inv_jacobian(int, int, Point);
};

/*! \class EFToNearPolePi
 * \brief Coordinate conversion from EF coordinates to stereographic near pole theta=pi
 */
class EFToNearPolePi : public CoordinateConversion
{
public:
	EFToNearPolePi();
	~EFToNearPolePi();
	
	Point convertPoint(Point);
	double jacobian(int, int, Point);
	double inv_jacobian(int, int, Point);
};

/*! \class NearPolePiToEF
 * \brief Coordinate conversion from stereographic coordinates near pole theta=pi to EF
 */
class NearPolePiToEF : public CoordinateConversion
{
public:
	NearPolePiToEF();
	~NearPolePiToEF();
	
	Point convertPoint(Point);
	double jacobian(int, int, Point);
	double inv_jacobian(int, int, Point);
};

/*! \class NearPoleToNearPole
 * \brief Coordinate conversion from stereographic coordinates near one pole to another
 */
class NearPoleToNearPole : public CoordinateConversion
{
public:
	NearPoleToNearPole();
	~NearPoleToNearPole();
	
	Point convertPoint(Point);
	double jacobian(int, int, Point);
	double inv_jacobian(int, int, Point);
};

#endif

