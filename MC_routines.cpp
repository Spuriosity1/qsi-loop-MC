#include "MC_routines.hpp"
#include "chain.hpp"
#include <cstring>
#include <random>

void initialise_spins(Lattice &lat, rng_t &gen) {
    // approximate half up half down initialisation
    std::bernoulli_distribution d(0.5);
    for (auto &s : lat.get_links()) {
        s.sz = 2*static_cast<int>(d(gen)) -1;
    }
}
double local_field(const Spin *s) {
    double h = 0.0;
    for (const auto &bond : s->nn_bonds) {
        if (bond.spin) {
            h += bond.Jzz * bond.spin->sz; // Jzz sz_i sz_j
        }
    }
    return h;
}

uint64_t apply_local_move(Spin *s, rng_t &gen, const MC_settings &settings) {

    double h = local_field(s);
    double dE = -2.0 * h * s->sz; // ΔE = E_after - E_before = 2 * h * sz

    if (settings.T <= 0) {
        if (dE > 0) {
            return 0; // reject
        }
    } else {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        if (dE > 0 && dist(gen) >= std::exp(-dE / settings.T)) {
            return 0; // reject
        }
    }
    s->sz *= -1; // accept
    return 1;
}

bool find_alternating_loop(Spin *start, std::vector<Spin *> &loop,
                                                      std::unordered_set<Spin *> &visited, rng_t &gen,
                                                      const MC_settings &settings) {

    // Step 1: pick a random direction from the up tetrahedron (sl = 0)
    std::vector<Spin *> candidates;

    for (int i = 0; i < 3; ++i) {
        const auto &bond = start->nn_bonds[i];
        Spin *neighbor = bond.spin;
#ifndef NDEBUG
        assert(neighbor);

        // Check that the two spins share a tetrahedron with sl = 0
        auto shared = start->boundary & neighbor->boundary;
        assert(shared.size() == 1);
        auto tetra = static_cast<Tetra *>(shared[0]);
        assert(tetra->sl == 0);
#endif
        if (antialigned(start, neighbor)) {
            candidates.push_back(neighbor);
        }
    }

    // If no valid neighbors found, no loop can be initiated from this spin
    if (candidates.empty())
        return false;

    std::uniform_int_distribution<size_t> choose(0, candidates.size() - 1);
    Spin *first = candidates[choose(gen)];

    assert(loop.empty());
    loop.push_back(start);
    loop.push_back(first);
    visited.insert(start);
    visited.insert(first);

    std::array<size_t, 3> idx = {0, 1, 2};

    // worm update
    sl_t tetra_sl = 0;
    while (loop.size() < settings.max_loop_len) {
        Spin *current = loop.back();
        tetra_sl ^= 1; // we change tetra sl every hop (faster than checking by hand)
        std::shuffle(idx.begin(), idx.end(), gen);
	
		int i=0;
        for (i = 0; i < 3; i++) {
            auto bond = current->nn_bonds[3 * tetra_sl + idx[i]];
            Spin *next = bond.spin;
#ifndef NDEBUG
            assert(next);

            auto shared = current->boundary & next->boundary;
            assert(shared.size() == 1);
            auto tetra = static_cast<Tetra *>(shared[0]);
            assert(tetra->sl == tetra_sl);

            if (loop.size() > 1) {
                assert(next != loop[loop.size() - 2]);
            } // bipartiteness should make this impossible
#endif
            if (visited.count(next))
                continue;
            if (!antialigned(current, next))
                continue;

            if (next == start && loop.size() >= 4) {
                // we did it
                return true;
            }

            loop.push_back(next);
            visited.insert(next);
			break;
        }
		if (i == 3) {
			// i the for loop reached the end, which should not happen in a 2I2O
			return false;
		}
    }

    return false;
}

double total_local_energy(const std::vector<Spin *> &loop) {
    // sum runs over all pairs

    auto tetra_0 = loop[0]->boundary & loop[1]->boundary;
    int tetra_sl = static_cast<Tetra *>(tetra_0[0])->sl;
    assert(tetra_sl == 0);

    double h = 0;

    for (size_t i = 0; i < loop.size(); i++) {
        auto s1 = loop[i];
        auto s2 = loop[(i + 1) % loop.size()];

        assert(antialigned(s1, s2));

        // find the non-[s1,s2] bonds
        //
        // nn_bonds[i] lies on sublattice (sl of this spin) + 1 + i % 4 joined by A
        // tetra (i.e. s->boundary[0]) nn_bonds[3+i] lies on sublattice (sl of this
        // spin) + 1 + i % 4 joined by B tetra (i.e. s->coboundary[0])
        /*
          *      s1
          *       .\ \
          *       . \   \
          *       .  s4   s3
          *       . /   /
          *       ./ /
          *      s2
          */

        for (int j = 0; j < 3; j++) {
            auto sn = s1->nn_bonds[3 * tetra_sl + j];
            if (sn.spin == s2) {
                continue;
            }
            h += sn.Jzz * sn.spin->sz;

            sn = s2->nn_bonds[3 * tetra_sl + j];
            assert(sn.spin != s1);
            h -= sn.Jzz * sn.spin->sz;
        }

        tetra_sl ^= 1;
        h *= -1;
    }
    return h * loop[0]->sz;
}

bool attempt_loop_move(Spin *root, rng_t &gen, const MC_settings &settings) {

    std::unordered_set<Spin *> visited;
    std::vector<Spin *> loop;

    if (!find_alternating_loop(root, loop, visited, gen, settings)) {
        return false; // no valid loop found
    }

    double E_before = total_local_energy(loop);

    // Flip the loop spins
    for (Spin *s : loop)
        s->sz *= -1;

    double E_after = total_local_energy(loop);
    double dE = E_after - E_before;

    bool accepted = false;

    if (settings.T <= 0) {
        accepted = dE <= 0;
    } else {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        accepted = (dE <= 0) || (dist(gen) < std::exp(-dE / settings.T));
    }

    if (!accepted) {
        // Revert flip
        for (Spin *s : loop)
            s->sz *= -1;
    }

    return accepted;
}

uint64_t mc_sweep_local(std::vector<Spin *> &spins, rng_t &gen,
                                                const MC_settings &settings) {
    uint64_t success = 0;
    if (settings.sweep_strategy == MC_settings::sweep::random) {
        std::uniform_int_distribution<size_t> pick(0, spins.size() - 1);
        for (size_t i = 0; i < spins.size(); ++i) {
            success += apply_local_move(spins[pick(gen)], gen, settings);
        }
    } else if (settings.sweep_strategy == MC_settings::sweep::sequential) {
        for (size_t i = 0; i < spins.size(); ++i) {
            success += apply_local_move(spins[i], gen, settings);
        }
    } else {
        throw std::logic_error("Unimplemented MC sweep strategy");
    }
    return success;
}

uint64_t mc_sweep_loop(std::vector<Spin *> &spins, rng_t &gen,
                                              const MC_settings &settings) {
    uint64_t success = 0;
    if (settings.sweep_strategy == MC_settings::sweep::random) {
        std::uniform_int_distribution<size_t> pick(0, spins.size() - 1);
        for (size_t i = 0; i < spins.size(); i++) {
            success += attempt_loop_move(spins[pick(gen)], gen, settings);
        }
    } else if (settings.sweep_strategy == MC_settings::sweep::sequential) {
        for (size_t i = 0; i < spins.size(); i++) {
            success += attempt_loop_move(spins[i], gen, settings);
        }
    } else {
        throw std::logic_error("Unimplemented MC sweep strategy");
    }
    return success;
}

