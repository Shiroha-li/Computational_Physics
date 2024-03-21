#ifndef IsingSystem_hpp
#define IsingSystem_hpp

#include<iostream>
#include<cassert>

class IsingSpin{
    private:
    // sz属性可以用于IsingSpin和SpinOnLattice的对象，表示该点处的spin state
    int sz;

    public:
    // initialize sz
    IsingSpin() {sz = 1;};
    ~IsingSpin() {};

    // Basic functions:display sz;set up/down/value;flip
    int _sz() const { return sz; };
    void set_up() {sz = 1;};
    void set_down() {sz = -1;};
    void set_sz(int sz_spec) {
        assert(sz_spec ==1 || sz_spec == -1);
        sz = sz_spec;
    };
    void flip() {sz *= -1;};
};

#endif