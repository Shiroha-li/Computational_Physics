#include<iostream>
#include<vector>
#include"spin.hpp"
#include"Ising_System.hpp"
using namespace std;
int main(int argc, const char * argv[]) 
{
    /*IsingSpin spin;
    cout << "Spin state (initiate): " << spin._sz() << endl;

    spin.flip();
    cout << "Spin state (initiate): " << spin._sz() << endl;

    spin.flip();
    cout << "Spin state (initiate): " << spin._sz() << endl;
    */
    cout << "How many particles?" << endl;
    int number = 0;
    cin >> number;
    const int n_spin = number;
    /*vector<IsingSpin> spin_array(n_spin);
    
    const double J = -1;
    double energy = 0;
    for (int i = 1;i<n_spin;i++) {
        energy += spin_array[i-1]._sz() * spin_array[i]._sz();
    }
    energy += spin_array[n_spin - 1]._sz() * spin_array[0]._sz();
    energy *= J;
    */

    vector<double> beta = {0,1};
    IsingSystem ISwithTenParticles(n_spin,beta);
    int Configuration[] = {7,77,777};
    for (int i = 0; i < 3; i++){
        ISwithTenParticles.set_state_by_code(Configuration[i]);
        cout << "When spin configuration is " << Configuration[i] << endl;
        cout<< "Magnetization = " << ISwithTenParticles.eval_mz() << endl;
        cout << "Energy = " << ISwithTenParticles.eval_energy_1D() << endl;
    }
    return 0;
}