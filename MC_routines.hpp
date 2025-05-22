#pragma once

#include "struct.hpp"
#include <cassert>
#include <cell_geometry.hpp>
#include <argparse.hpp>
#include <chain.hpp>
#include <cstddef>
#include <cstdint>
#include <lattice_IO.hpp>
#include <preset_cellspecs.hpp>
#include <UnitCellSpecifier.hpp>
#include <fftw3.h>
#include <string>
#include <XoshiroCpp.hpp>
#include <unordered_set>


// the RNG engine of choice
// typedef std::mt19937 rng_t;
typedef XoshiroCpp::Xoshiro256PlusPlus rng_t;

void initialise_spins(Lattice &lat, rng_t &gen);

struct MC_settings {
    enum class sweep { sequential, random };
    sweep sweep_strategy;
    double T=0;
    size_t max_loop_len = 1<<31; // size cutoff for the loop-find algo
    uint64_t n_sweeps; // total number of sweeps to do
    uint64_t sweep_n_local; // number of local updates per sweep
    uint64_t sweep_n_loop; // number of local updates per sweep                           
};



// Returns sum(bond.Jzz * s_i) for s_i in nn
double local_field(const Spin *s);

inline double local_energy(const Spin* s) {
    return local_field(s) * s->sz;
}

uint64_t apply_local_move(Spin *s, rng_t &gen, const MC_settings &settings);

inline bool antialigned(const Spin* a, const Spin* b) {
    return a->sz * b->sz < 0;
}

/** Attempts to find a closed, flippable loop centred on 'start'.
 * Uses DFS algorithm as follows:
 * 1. 
 *
 */
bool find_alternating_loop(Spin *start, std::vector<Spin *> &loop,
                           std::unordered_set<Spin *> &visited, rng_t &gen,
                           const MC_settings &settings);

double total_local_energy(const std::vector<Spin *> &loop);

// Finds a random closed loop of alternating spins and attempts to flip it.
bool attempt_loop_move(Spin *root, rng_t &gen, const MC_settings &settings);

uint64_t mc_sweep_local(std::vector<Spin *> &spins, rng_t &gen,
                        const MC_settings &settings);

uint64_t mc_sweep_loop(std::vector<Spin *> &spins, rng_t &gen,
                       const MC_settings &settings);




