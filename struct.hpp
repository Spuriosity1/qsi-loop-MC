#pragma once
#include <cassert>
#include <cell_geometry.hpp>
#include <argparse.hpp>
#include <chain.hpp>
#include <cstddef>
#include <lattice_IO.hpp>
#include <preset_cellspecs.hpp>
#include <UnitCellSpecifier.hpp>
#include <fftw3.h>
#include <XoshiroCpp.hpp>


using namespace CellGeometry;
using namespace nlohmann;

extern "C" {
#include <hdf5.h>
}

struct Spin; // fwd declaration

struct Tetra : public Cell<0> {
    // the NN couplings, understood as a symmetric Jzz[i, j].
    // double Jzz[4*4] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int sl;
};


struct Bond {
    double Jzz; // the Ising coupling
    Spin* spin;
};

struct Spin : public Cell<1> {
    int sz; // Ising spin +-1
    Bond nn_bonds[6];
    // For i=0,1,2; nn_bonds[i] is over the 'A' tetrahedron, 'B' for 3,4,5.
    // if "sl2" is the sl of a nn, it can be accessed at
    // nn_bonds[tetra_sl*3 + neighbor_index[sl1][sl2]] is guaranteed to lie on sublattice
    // neighbor_index[lat.sl_of_link(*this)]
    constexpr static const int neighbor_index[4][4] = {
        { 128, 0, 1, 2 },
        { 0, 128, 1, 2 },
        { 0, 1, 128, 2 },
        { 0, 1, 2, 128 }
    };  
};

struct Plaq : public Cell<2> {
};

struct Vol : public Cell<3> {
};


template <int order>
inline std::vector<Cell<order>*> operator&(const Chain<order>& a, const Chain<order>& b){
    std::vector<Cell<order>*> retval;
    for (auto [s,_] : a){
        if (b.count(s)) {
            retval.push_back(s);
        }
    }
    return retval;
} 


typedef PeriodicPlaqLattice<Tetra, Spin, Plaq> Lattice;

////////////////////////////////////////////////////////////////////////////////
//////// REGISTRATION OF OBJECTS AND CACHING ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


template <typename func>
void setup_couplings(Lattice& lat, func Jzz_func){
    // set up the NN shortcuts

    for (auto& [_, t] : lat.points){
        t->sl = lat.primitive_spec.sl_of_point(t->position);
        for (auto [l1, m] : t->coboundary){ 
            for (auto [l2, m] : t->coboundary){
                if (l1 <= l2) continue; // avoids double counting
                auto s1 = static_cast<Spin*>(l1);
                auto s2 = static_cast<Spin*>(l2);

                auto s1_sl= lat.primitive_spec.sl_of_link(s1->position);
                auto s2_sl= lat.primitive_spec.sl_of_link(s2->position);

                assert(s1_sl != s2_sl); 
                auto Jzz = Jzz_func(s1_sl,s2_sl);

                // t->Jzz[4*s1_sl + s2_sl] = Jzz; 
                // t->Jzz[4*s2_sl + s1_sl] = Jzz; 
                s1->nn_bonds[3*t->sl + Spin::neighbor_index[s1_sl][s2_sl]] = Bond({Jzz, s2});
                s2->nn_bonds[3*t->sl + Spin::neighbor_index[s2_sl][s1_sl]] = Bond({Jzz, s1});

            }
        }
    } 
}


