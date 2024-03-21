#include <iostream>
#include <vector>
#include <cmath>
#include "../Ising_System.hpp"
using namespace std;

class IsingSystem_Square : public IsingSystem {
    private:
        vector<int> system_size;
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
        IsingSystem_Square(const vector<int> system_size_spec) :
            IsingSystem(system_size_spec[0] * system_size_spec[1]),
            system_size(system_size_spec) { setup_NN(); }; // initialize NN
        ~IsingSystem_Square() {};
        // Encode (这里只考虑了系统是二维的情况)
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
        // set_state接收表示spin state的bool数组，并给相应对象的sz属性赋值
        void set_state(vector<bool> state) { 
            if (state.size() == n_spins) {
                for (int i = 0; i < n_spins ; ++i) {
                    if (state[i] == 1) 
                        set_up_spin(i);
                    else
                        set_dw_spin(i);
                }
            }
            else 
                cout << "输入的state数组大小与n_spins不符" ;
        };
        double eval_Mz() { return eval_mz(); }; // 这行可以不要
        double eval_energy() {
            double Energy = 0;
            for (int site_idx = 0; site_idx < n_spins; site_idx++ ) {
                for (int bond_idx = 0; bond_idx < 4; bond_idx++ ) {
                    Energy += _sz(site_idx) * _sz(NN(site_idx, bond_idx));
                }
            }
            Energy *= J/2;
            return Energy;
        };
};
