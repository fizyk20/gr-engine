#ifndef __KERR_H__
#define __KERR_H__

/*! \file kerr.h
 * \brief Implementation of the Kerr metric
 */

#include "geometry.h"
#include "kerr_coords.h"

/*! \class KerrManifold
 * \brief Class representing a Kerr manifold
 */
class KerrManifold : public Manifold
{
	double M, a;
public:
	KerrManifold(double, double);
	~KerrManifold();
	
	double getMass();
	double getAngMomentum();
	void setMass(double _M);
	void setAngMomentum(double _a);
	
	int recommendCoordSystem(Point);
};

/*! \class KerrEFMetric
 * \brief The Schwarzschild metric in Eddington-Finkelstein coordinates
 */
class KerrEFMetric : public Metric
{
	KerrManifold* m;
protected:
	double _g(int, int, Point);
	double _invg(int, int, Point);
	double _christoffel(int, int, int, Point);
	
public:
	enum { coordU = 0, coordR = 1, coordTheta = 2, coordPhi = 3 };
	
	KerrEFMetric(int cS, KerrManifold* _m);
	~KerrEFMetric();
};

/*! \class KerrNearPoleMetric
 * \brief The Schwarzschild metric valid near spherical poles - in "stereographic" coordinates
 */
class KerrNearPoleMetric : public Metric
{
	KerrManifold* m;
protected:
	double _g(int, int, Point);
	double _invg(int, int, Point);
	double _christoffel(int, int, int, Point);
	
public:
	enum { coordU = 0, coordR = 1, coordX = 2, coordY = 3 };
	
	KerrNearPoleMetric(int cS, KerrManifold* _m);
	~KerrNearPoleMetric();
};

#endif

