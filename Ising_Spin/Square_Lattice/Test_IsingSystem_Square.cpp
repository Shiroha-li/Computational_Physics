#include<catch2/catch_test_macros.hpp>
#include<catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include <vector>
//#include "SpinOnLattice.hpp" //Test Case 1
//#include "../Ising_System.hpp" //Test Case 2
#include "IsingSystem_Square.hpp" //Test Case 3
using namespace std;

TEST_CASE("IsingSpinOnLattice","[single spin]"){
    IsingSpinOnLattice spin;
    SECTION("spin position/sublattice (initial)"){
        constexpr int dim=2;
        spin.set_dim(dim);
        REQUIRE(spin._position() == vector<int>({0,0}));
        REQUIRE(spin._NN() == vector<int>({-1}));
        REQUIRE(spin._NN(0) == -1);
    }
};

TEST_CASE("IsingSystem", "[examples of 10 spins]") {
    // 实例化了一个有10个粒子的对象model
    const int n_spin = 10;
    IsingSystem model(n_spin);

    SECTION("basics"){}

    SECTION("spin state code #7"){}

    SECTION("spin state code #77"){}

    SECTION("spin state code #777"){}

    SECTION("spin position initialization"){
        constexpr int dim = 2;
        model.set_dim(dim);
        for (int i = 0 ; i < n_spin ; i++) {
            REQUIRE(model._spin_position(i) == vector<int>({0,0}));
            REQUIRE(model._spin_NN(i) == vector<int>({-1}));
            REQUIRE(model._spin_NN(i, 0) == -1);
        }
    }
};

TEST_CASE("IsingSystem_Square", "[examples of 6 x 6 spins]") {
    const vector<int> system_size= {6, 6};
    IsingSystem_Square model(system_size);
    
    SECTION("basics") {
        REQUIRE(model._n_spins()== 36);
        REQUIRE_THAT(model._J(), Catch::Matchers::WithinULP(-1.0, 4));
    }

    SECTION("site index") {
        const vector<int> lattice_coordinate ={3, 4};
        REQUIRE(model.site_index(lattice_coordinate) == 27); // 4*6+3=27
        REQUIRE(model.lattice_coordinate(27) == lattice_coordinate); 
    }

    SECTION("neighboring site coordinates") {
        const vector<int> lattice_coordinate ={3, 3};
        const vector<int> lattice_coordinate_pos_x={4, 3};
        const vector<int> lattice_coordinate_pos_y={3, 4};
        const vector<int> lattice_coordinate_neg_x={2, 3};
        const vector<int> lattice_coordinate_neg_y={3, 2};
        REQUIRE(model.shift_pos_x(lattice_coordinate) == lattice_coordinate_pos_x);
        REQUIRE(model.shift_pos_y(lattice_coordinate) == lattice_coordinate_pos_y);
        REQUIRE(model.shift_neg_x(lattice_coordinate) == lattice_coordinate_neg_x);
        REQUIRE(model.shift_neg_y(lattice_coordinate) == lattice_coordinate_neg_y);
    }

    SECTION("connectivity") {
        constexpr int i = 21;
        REQUIRE(model.NN(i, 0) == 22);
        REQUIRE(model.NN(i, 1) == 27);
        REQUIRE(model.NN(i, 2) == 20);
        REQUIRE(model.NN(i, 3) == 15);

        constexpr int j = 0;
        REQUIRE(model.NN(j, 0) == 1);
        REQUIRE(model.NN(j, 1) == 6);
        REQUIRE(model.NN(j, 2) == 5);
        REQUIRE(model.NN(j, 3) == 30);
    }

    SECTION("'pi' state : magnetization and energy") {
        vector<bool>
            state({1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0});
        model.set_state(state);
        REQUIRE(model.eval_Mz() == 2);
        REQUIRE_THAT(model.eval_energy(), Catch::Matchers::WithinULP(-4.0, 4));
    }
};