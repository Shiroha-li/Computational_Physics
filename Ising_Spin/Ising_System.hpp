#include<iostream>
#include<vector>
#include<cmath>
#include"Square_Lattice/SpinOnLattice.hpp"
using namespace std;

class IsingSystem {
protected:
    const double J;
    const int n_spins;
    const long long maxrep_state;
    vector<IsingSpinOnLattice> spin;

public:
    IsingSystem(const int n_spins_spec) : J(-1.0), n_spins(n_spins_spec), maxrep_state(static_cast<long long>(pow(2,n_spins))-1) {
        spin.resize(n_spins);
    };
    virtual ~IsingSystem() {};

    double _J() const { return J;};
    int _n_spins() const { return n_spins;};
    long long _maxrep_state() const { return maxrep_state;};

    int _sz(const int site_idx) const { return spin[site_idx]._sz();}
    void set_up_spin(const int site_idx){spin[site_idx].set_up();};
    void set_dw_spin(const int site_idx){spin[site_idx].set_down();};
    void set_spin(const int site_idx, int s_spec)
        { spin[site_idx].set_sz(s_spec);};
    void flip_spin(const int site_idx){ spin[site_idx].flip();};
    void set_state_by_code(long long rep_state){
        for (int i = 0 ; i < n_spins ; i++){
            if (rep_state & 1)
                set_up_spin(i);
            else
                set_dw_spin(i);
            rep_state >>= 1;
        }
    };
    double eval_mz()const{
        double Magnetization=0;
        for(int t = 0 ; t < n_spins ; t++){
            Magnetization += _sz(t);
        }
        return Magnetization;
    };
    double eval_energy_1D()const{
        double H=0;
        if( n_spins ==1 &&  _sz(0) ) 
            H=0.5;
        else{
            for(int i=0;i<n_spins-1;i++){
                H += _sz(i)*_sz(i+1);           
            };
            H += _sz(n_spins-1) * _sz(0);
            H *= J;
        }
        return H;
    };
  
    // Additional part in Lattice case
    void set_dim(int dim) { for (auto& each: spin) each.set_dim(dim); };
    vector<int> _spin_position(const int site_idx) const {
        return spin[site_idx]._position();
    };
    vector<int> _spin_NN(const int site_idx) const {
        return spin[site_idx]._NN();
    };
    int _spin_NN(const int site_idx, const int bond_idx) const { 
        return spin[site_idx]._NN(bond_idx);
    };
};
