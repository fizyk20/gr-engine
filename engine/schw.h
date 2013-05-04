#ifndef __SCHW_H__
#define __SCHW_H__

/*! \file schw.h
 * \brief Implementation of the Schwarzschild metric
 */

#include "geometry.h"
#include "kerr_coords.h"

/*! \class SchwManifold
 * \brief Class representing a Schwarzschild manifold
 */
class SchwManifold : public Manifold
{
	double M;
public:
	SchwManifold(double);
	~SchwManifold();
	
	double getMass();
	void setMass(double _M);
	
	int recommendCoordSystem(Point);
};

/*! \class SchwEFMetric
 * \brief The Schwarzschild metric in Eddington-Finkelstein coordinates
 */
class SchwEFMetric : public Metric
{
	SchwManifold* m;
protected:
	double _g(int, int, Point);
	double _invg(int, int, Point);
	double _christoffel(int, int, int, Point);
	
public:
	enum { coordU = 0, coordR = 1, coordTheta = 2, coordPhi = 3 };
	
	SchwEFMetric(int cS, SchwManifold* _m);
	~SchwEFMetric();
};

/*! \class SchwNearPoleMetric
 * \brief The Schwarzschild metric valid near spherical poles - in "stereographic" coordinates
 */
class SchwNearPoleMetric : public Metric
{
	SchwManifold* m;
protected:
	double _g(int, int, Point);
	double _invg(int, int, Point);
	double _christoffel(int, int, int, Point);
	
public:
	enum { coordU = 0, coordR = 1, coordX = 2, coordY = 3 };
	
	SchwNearPoleMetric(int cS, SchwManifold* _m);
	~SchwNearPoleMetric();
};

#endif

