#pragma once
#include "rationalmath.hpp"
#include "struct.hpp"


////////////////////////////////////////////////////////////////////////////////
///////   OBSERVABLES                   ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T>
T convert_rational(const rational::Rational& r){
    return static_cast<T>(r.num)/r.denom;
}


// class responsible for calculating \sum_i \sigma_i e^{iq \cdot r}
class SSF_manager {
    const Lattice& lat;

    fftw_plan plans[4]; 
    // temporary storage
    fftw_complex* fft_out;
    double* spin_field[4]; 

    // the accumulator for <Sz(q)Sz(-q)> (note it is real)
    double* S_q_acc[4]; // indexed by sublattice
    uint64_t n_samples=0;

    hsize_t S_q_acc_dims[2];

public:

    SSF_manager(const Lattice& _lat);
    ~SSF_manager();

    void store();
    void save_to_hdf5(const std::string& filename);
};
