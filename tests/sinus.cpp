#include "../numeric.h"
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
    RK4Integrator integrator;
    
    StateVector state0;
    state0.push_back(0.0);
    state0.push_back(1.0);
    
    integrator.setStepSize(0.04);
    integrator.setState(state0);
    integrator.setDiffEquation(&eq);
    
    int i;
    ofstream fout("test.dat");
    for(i=0; i<400; i++)
    {
    	StateVector v = integrator.getState();
    	cout << "x = " << v[0] << "       v = " << v[1] << endl;
    	fout << integrator.getT() << "\t" << v[0] << "\t" << v[1] << endl;
    	integrator.next();
    }
    fout.close();
    
    return 0;
}
