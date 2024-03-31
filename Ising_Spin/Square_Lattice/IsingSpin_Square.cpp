#include <iostream>
#include "IsingSystem_Square.hpp"

using namespace std;

int main()
{
    double start = 0.05; 
    double end = 4.0; 
    double step = 0.05; 

    vector<double> Temperature;
   
    for (double i = start; i < end; i += step) {
        Temperature.push_back(i); 
    }
    vector<double> beta(Temperature.size());
    
    for (int i = 0; i < Temperature.size(); i++) {
        beta[i] = 1.0 / Temperature[i];
    }

    for (int i = 2; i < 5; i++){
        vector<int> system_size = { i, i };
        IsingSystem_Square model(system_size, beta);
        model.exact();
        model.print_exact();
    }
    // model.weight_unnormalized(1,0);

    return 0;

};