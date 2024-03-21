#include <iostream>
#include <vector>
#include <cmath>
#include "../spin.hpp"
using namespace std;

//First Inheritance
class IsingSpinOnLattice : public IsingSpin {
    private:
        vector<int> position;
        vector<int> NN;
    public:
        // Constructor
        IsingSpinOnLattice() : NN({-1}) {};
        // Destructor
        ~IsingSpinOnLattice(){};

        // Set value , 将dim个0赋值给position数组
        void set_dim(int dim){
            position.assign(dim,0);
        };

        //Display the value of position and NN
        vector<int> _position() const { return position; };
        vector<int> _NN() const { return NN; };
        int _NN(const int bond_idx) const { return NN[bond_idx]; }; 

        void set_NN(const int bond_idx, const int site_idx) { NN[bond_idx] = site_idx; };

};