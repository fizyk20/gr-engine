#include "kerr_coords.h"
#include <math.h>

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
	if(p.getCoordSystem() != EF) throw "EFToNearPole0: invalid coordinate system.";
	Point result(NearPole0);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = tan(p[2]/2)*cos(p[3]);
	result[3] = tan(p[2]/2)*sin(p[3]);
	
	return result;
}

double EFToNearPole0::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != EF) throw "EFToNearPole0: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != EF) throw "EFToNearPole0: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != EF) throw "EFToNearPolePi: invalid coordinate system.";
	Point result(NearPolePi);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = tan((M_PI-p[2])/2)*cos(p[3]);
	result[3] = tan((M_PI-p[2])/2)*sin(p[3]);
	
	return result;
}

double EFToNearPolePi::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != EF) throw "EFToNearPolePi: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != EF) throw "EFToNearPolePi: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	Point result(EF);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = 2*atan(sqrt(p[2]*p[2]+p[3]*p[3]));
	result[3] = atan2(p[3],p[2]);
	
	return result;
}

double NearPole0ToEF::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != NearPole0) throw "NearPole0ToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	Point result(EF);
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = M_PI - 2*atan(sqrt(p[2]*p[2]+p[3]*p[3]));
	result[3] = atan2(p[3],p[2]);
	
	return result;
}

double NearPolePiToEF::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != NearPolePi) throw "NearPolePiToEF: invalid coordinate system.";
	
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
	if(p.getCoordSystem() != NearPole0 && p.getCoordSystem() != NearPolePi) 
		throw "NearPole0ToNearPolePi: invalid coordinate system.";
	
	Point result(p.getCoordSystem() == NearPole0 ? NearPolePi : NearPole0);
	
	double l = p[2]*p[2] + p[3]*p[3];
	
	result[0] = p[0];
	result[1] = p[1];
	result[2] = p[2]/l;
	result[3] = p[3]/l;
	
	return result;
}

double NearPoleToNearPole::jacobian(int i, int j, Point p)
{
	if(p.getCoordSystem() != NearPole0 && p.getCoordSystem() != NearPolePi) 
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
	if(p.getCoordSystem() != NearPole0 && p.getCoordSystem() != NearPolePi) 
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

