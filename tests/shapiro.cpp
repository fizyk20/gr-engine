#include "../engine/particle.h"
#include "../engine/dpintegrator.h"
#include "../engine/schw.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
using namespace std;

double u(double t, double r, double M)
{
	return t + r + 2*M*log(0.5*(r-2*M)/M);
}

double t(double u, double r, double M)
{
	return u - r - 2*M*log(0.5*(r-2*M)/M);
}

int main(int argc, char** argv)
{
	double M = 4.9e-6;	//Słońce
	double d = 2.33;
	double yE = 498.67;
	double yV = 370.7;
	
	cout << "The program calculates the Shapiro delay for a ray of light." << endl;
	cout << "Usage: shapiro [d [yV [yE]]]" << endl;
	cout << "d - the perihelion of the ray in light-seconds" << endl;
	cout << "yV - the y coordinate of Venus (= sqrt(rV^2-d^2), where rV is the distance of Venus from the Sun) in light-seconds" << endl;
	cout << "yE - the y coordinate of Earth (= sqrt(rE^2-d^2), where rE is the distance of Earth from the Sun) in light-seconds" << endl;
	cout << "Defaults: d = 2.33, yV = 370.7, yE = 498.67" << endl << endl;
	
	switch(argc)
	{
	case 1:
		cout << "Data set to default." << endl;
		break;
	case 2:
		d = atof(argv[1]);
		break;
	case 3:
		d = atof(argv[1]);
		yV = atof(argv[2]);
		break;
	case 4:
		d = atof(argv[1]);
		yV = atof(argv[2]);
		yE = atof(argv[3]);
		break;
	}
	cout << "d = " << d << endl << "yE = " << yE << endl << "yV = " << yV << endl;
	double u0 = sqrt(d*d*d/(d-2*M));
	double rE = sqrt(d*d + yE*yE);
	double rV = sqrt(d*d + yV*yV);
	
	double t_flat = 2*(yE + yV);
	
	SchwManifold schw(M);
	vector4 u_init = vector4(u0, 0.0, 0.0, 1.0);
	
	Particle photon1(&schw, Point(EF, u(0.0, d, M), d, M_PI/2, 0.0), u_init);
	Particle photon2(&schw, Point(EF, u(0.0, d, M), d, M_PI/2, 0.0), vector4()-u_init);
	
	DPIntegrator dp(1e-12);
	photon1.setIntegrator(&dp);
	photon2.setIntegrator(&dp);
	
	double t1, t2;
	vector4 lastPos, pos;
	
	cout << "Propagation of photon 1..." << endl;
	
	while(photon1.getPos()[1] < rE)
	{
		lastPos = photon1.getPos().toVector4();
		photon1.propagate();
	}
	pos = photon1.getPos().toVector4();
	cout << "r = " << pos[1] << endl;
	
	lastPos += (pos - lastPos)/(pos[1] - lastPos[1])*(rE - lastPos[1]);
	t1 = t(lastPos[0], lastPos[1], M);
	
	dp.resetStepSize();
	
	cout << "Propagation of photon 2..." << endl;
	
	while(photon2.getPos()[1] < rV)
	{
		lastPos = photon2.getPos().toVector4();
		photon2.propagate();
	}
	pos = photon2.getPos().toVector4();
	
	lastPos += (pos - lastPos)/(pos[1] - lastPos[1])*(rV - lastPos[1]);
	t2 = t(lastPos[0], lastPos[1], M);
	
	cout << "Propagation finished." << endl;
	cout.precision(8);
	cout << "t1 = " << t1 << endl;
	cout << "t2 = " << t2 << endl;
	
	double dt = (t1 - t2)*2;
	dt *= sqrt(1.0 - 2*M/rE);	//time dilation from the Sun
	cout << "dt = " << dt << endl;
	cout << "delay = " << dt - t_flat << endl;
	
	return 0;
}
