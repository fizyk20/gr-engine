#ifndef __KERR_COORDS_H__
#define __KERR_COORDS_H__

/*! \file kerr_coords.h
 * \brief Coordinate conversions for Kerr / Schwarzschild metrics
 */

#include "geometry.h"

#define EF 0			/// Eddington-Finkelstein coordinates
#define NearPole0 1		/// Near-pole ("stereographic") coordinates near theta=0
#define NearPolePi 2	/// Near-pole ("stereographic") coordinates near theta=pi

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

