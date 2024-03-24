#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "../Ising_System.hpp"
using namespace std;

class IsingSystem_Square : public IsingSystem {
    private:
        vector<int> system_size;
        vector<double> internal_E;
        vector<double> internal_E_sq;
        vector<double> M_sq;
        vector<double> Z;
        vector<double> C;
        double ground_state_energy;

        // SpinOnLattice的对象spin[site_idx]调用SpinOnLattice的方法set_NN
        // set_NN接受bond_idx和对应的site_idx作为参数，将NN[bond_idx]赋值为site_index(shift...)
        // 经过构造函数初始化调用后，每一个spin[site_idx]对象的NN属性都成为一个含有4个相邻site_idx的数组
        // 并且SpinOnLattice的对象spin[site_idx]可以调用_NN方法来返回指定bond_idx处的site_idx
        void setup_NN() {
            for (int site_idx = 0 ; site_idx < n_spins ; site_idx++) {
                vector<int> r = lattice_coordinate(site_idx);
                spin[site_idx].set_NN(0, site_index(shift_pos_x(r)));
                spin[site_idx].set_NN(1, site_index(shift_pos_y(r)));
                spin[site_idx].set_NN(2, site_index(shift_neg_x(r)));
                spin[site_idx].set_NN(3, site_index(shift_neg_y(r)));
            }
        };

    public:
        // Constructor, take system_size_spec as input, then initialize system_size and calculate n_spins_spec for IsingSystem
        IsingSystem_Square(const vector<int> system_size_spec, const vector<double> beta_spec) :
            IsingSystem(system_size_spec[0] * system_size_spec[1], beta_spec),
            system_size(system_size_spec) {
                setup_NN(); // initialize NN
                internal_E.assign(beta.size(),0);
                internal_E_sq.assign(beta.size(),0);
                M_sq.assign(beta.size(),0);
                Z.assign(beta.size(),0);
                C.assign(beta.size(),0);
                ground_state_energy = eval_energy();
            }; 

        ~IsingSystem_Square() {};
        // Encode (这里只针对系统是二维矩形的情况)
        int site_index(const vector<int> lattice_coordinate) const { return lattice_coordinate[1] * system_size[0] + lattice_coordinate[0]; };
        // Decode
        vector<int> lattice_coordinate(int site_index) const { return vector<int>({site_index % system_size[1], site_index / system_size[1]}); };

        // Shifting to neighbouring site , r is for not changing the value of r_spec , % is for not exceeding the system_size
        vector<int> shift_pos_x(const vector<int> r_spec) const {
            vector<int> r(r_spec);
            r[0] = (r[0] + 1) % system_size[0];
            return r;
        };
        vector<int> shift_pos_y(const vector<int> r_spec) const{
            vector<int> r(r_spec);
            r[1] = (r[1] + 1) % system_size[1];
            return r;
        };
        vector<int> shift_neg_x(const vector<int> r_spec) const{
            vector<int> r(r_spec);
            r[0] = (r[0] - 1 + system_size[0]) % system_size[0];
            return r;
        };
        vector<int> shift_neg_y(const vector<int> r_spec) const{
            vector<int> r(r_spec);
            r[1] = (r[1] - 1 + system_size[1]) % system_size[1];
            return r;
        };

        // Getting the site index of a neighboring site
        /*
        int NN(const int site_idx, const int bond_idx) const {
            vector<int> r = lattice_coordinate(site_idx);
            switch(bond_idx){
                case 0:
                    return site_index(shift_pos_x(r));
                    break;
                case 1:
                    return site_index(shift_pos_y(r));
                    break;
                case 2:
                    return site_index(shift_neg_x(r));
                    break;
                case 3:
                    return site_index(shift_neg_y(r));
                    break;
                default:
                    return -1;
                }
        };
        */

        // Refactor，NN接受site_idx和bond_idx，返回相邻处的site_idx
        int NN(const int site_idx, const int bond_idx) const { return spin[site_idx]._NN(bond_idx); };

        // Evaluate M & E
        // 计算Magnetization可以直接用IsingSystem中的eval_mz()方法
        double eval_energy() const {
            double Energy = 0;
            for (int site_idx = 0; site_idx < n_spins; site_idx++ ) {
                for (int bond_idx = 0; bond_idx < 4; bond_idx++ ) {
                    Energy += _sz(site_idx) * _sz(NN(site_idx, bond_idx));
                }
            }
            Energy *= J/2;
            return Energy;
        };

        double ground_state() const {
            return ground_state_energy;
        };

        double weight_unnormalized(const size_t beta_idx) const {
            return exp(-beta[beta_idx]* ( eval_energy() - ground_state()));
        };

        double _exact_energy_Z(const size_t beta_idx) const {
            return weight_unnormalized(beta_idx);
        };

        double _exact_energy_q(const size_t beta_idx) const {
            return _exact_energy_Z(beta_idx) * eval_energy();
        }; 
    
        double _exact_energy_q_sq(const size_t beta_idx) const {
            return _exact_energy_q(beta_idx) * eval_energy();
        }; 
   
        double _exact_magz_Z(const size_t beta_idx) const {
            return weight_unnormalized(beta_idx);
        }; 

        double _exact_magz_q_sq(const size_t beta_idx) const {
            return _exact_magz_Z(beta_idx) * eval_mz() * eval_mz();
        }; 
    
        void exactly_evaluate_given() {
            for(int i = 0; i < beta.size(); i++) {
                internal_E[i] += _exact_energy_q(i);
                internal_E_sq[i] += _exact_energy_q_sq(i);
                Z[i] += weight_unnormalized(i);
                M_sq[i] += _exact_magz_q_sq(i);
            }
        };

        // For state in vector form
        void exactly_evaluate(const vector<bool>& state) {
            set_state(state),
            exactly_evaluate_given();
        };
    
        // For state in integer form
        void exactly_evaluate(const long long& rep_state) {
            std::vector<bool> state = state_by_code(rep_state);
            exactly_evaluate(state);
        };
    
        //going through all the state
        void exact() {
            long long rep_state = 0;
            while (rep_state <= maxrep_state) {
                exactly_evaluate(rep_state++);
            }
            normalize_direct();
        };

        void normalize_direct() {
            for(int i = 0; i<beta.size(); i++) {
                internal_E[i] *= 1/Z[i];
                internal_E_sq[i] *= 1/Z[i];
                M_sq[i] *= 1/Z[i];
                M_sq[i] /= n_spins * n_spins;
            }
            for(int i = 0; i<beta.size(); i++) {
                C[i] = beta[i] * beta[i] * (internal_E_sq[i]-internal_E[i]*internal_E[i]);
                C[i] /= n_spins;
            }
        };

        void print_exact() const {
            cout << "Specific Heat: ";
            for (double value : C) {
                cout << value << "  ";
            }
            cout << "end" << endl;
        
            cout << "Magnetization (Squared): ";
            for (double value : M_sq) {
                cout << value << "  ";
            }
            cout << "end" << endl;
        };

};
