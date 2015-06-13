#include "../engine/dpintegrator.h"
#include <iostream>
#include <fstream>
using namespace std;

class sinus_test : public DiffEq
{
    StateVector derivative(StateVector v)
    {
        if(v.size() != 3) throw StateLengthError();
        StateVector result;
        result.push_back(1.0); // dt/dt
        result.push_back(v[2]); // dx/dt
        result.push_back(-v[1]); // dv/dt
        
        return result;
    }
};

int main()
{
    sinus_test eq;
    DPIntegrator integrator(0.001);
    
    StateVector v;
    v.push_back(0.0); //t
    v.push_back(0.0); //x
    v.push_back(1.0); //v
    
    integrator.setStepSize(0.04);
    
    int i;
    ofstream fout("test_dp.dat");
    for(i=0; i<400; i++)
    {
    	cout << "step = " << integrator.getStepSize() << "   x = " << v[1] << "   v = " << v[2] << "   x^2 + v^2 = " << v[1]*v[1] + v[2]*v[2] << endl;
    	fout << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << endl;
    	v = integrator.next(v, &eq);
    }
    fout.close();
    
    return 0;
}
