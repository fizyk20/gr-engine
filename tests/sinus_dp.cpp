#include "../dpintegrator.h"
#include <iostream>
#include <fstream>
using namespace std;

class sinus_test : public DiffEq
{
    StateVector derivative(double t, StateVector v)
    {
        if(v.size() != 2) throw StateLengthError();
        StateVector result;
        result.push_back(v[1]);
        result.push_back(-v[0]);
        
        return result;
    }
};

int main()
{
    sinus_test eq;
    DPIntegrator integrator(0.001);
    
    StateVector state0;
    state0.push_back(0.0);
    state0.push_back(1.0);
    
    integrator.setStepSize(0.04);
    integrator.setState(state0);
    
    int i;
    ofstream fout("test_dp.dat");
    for(i=0; i<400; i++)
    {
    	StateVector v = integrator.getState();
    	cout << "step = " << integrator.getStepSize() << "   x = " << v[0] << "   v = " << v[1] << "   x^2 + v^2 = " << v[0]*v[0] + v[1]*v[1] << endl;
    	fout << integrator.getT() << "\t" << v[0] << "\t" << v[1] << "\t" << endl;
    	integrator.next(&eq);
    }
    fout.close();
    
    return 0;
}
