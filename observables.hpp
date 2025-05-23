#pragma once
#include "rationalmath.hpp"
#include "struct.hpp"
#include <cassert>

////////////////////////////////////////////////////////////////////////////////
///////   OBSERVABLES                   ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T> T convert_rational(const rational::Rational &r) {
    return static_cast<T>(r.num) / r.denom;
}

// abstract base class implementing the basic interface
class observable_manager {
    protected:
        const Lattice &lat;

        uint64_t n_samples = 0;

        hid_t open_hdf5(const std::string &filename);

        void write_data(hid_t file_id);

    public:
        observable_manager(const Lattice &_lat) : lat(_lat) { n_samples = 0; }

        void store() { n_samples++; }; // stores the current state of 'lat'
        void save_to_hdf5(const std::string &filename) {
            hid_t file_id = open_hdf5(filename);
            write_data(file_id);
            H5Fclose(file_id);
        }
};

// class responsible for calculating \sum_i \sigma_i e^{iq \cdot r}
class ssf_manager : public observable_manager {
        fftw_plan plans[4];
        // temporary storage
        fftw_complex *fft_out;
        double *spin_field[4];

        // the accumulator for <Sz(q)Sz(-q)> (note it is real)
        double *S_q_acc; // indexed by sublattice

        hsize_t S_q_acc_dims[3];

        void store_ssf();

    protected:
        void write_data(hid_t file_id);

    public:
        ssf_manager(const Lattice& _lat);
        ~ssf_manager();

        void store() {
            observable_manager::store();
            store_ssf();
        }

        void save_to_hdf5(const std::string &filename) {
            hid_t file_id = open_hdf5(filename);
            observable_manager::write_data(file_id);
            write_data(file_id);
            H5Fclose(file_id);
        }
};


// adds on raw <Sz> measurements
class Sz_manager : public ssf_manager {
    double *Sz_acc;
    void store_Sz(){
        for (auto& [J, l] : lat.links){
            auto s = static_cast<Spin*>(l);
            Sz_acc[J] += s->sz;
        }
    }

    protected:
    void write_data(hid_t file_id);

    public:
    Sz_manager(const Lattice& _lat) : ssf_manager(_lat) {
        // Sz_acc = new double[lat.num_primitive * 4];
    }
    ~Sz_manager(){
        // delete[] Sz_acc;
    }

    void store(){
        ssf_manager::store();
        // store_Sz();
    }

    void save_to_hdf5(const std::string &filename) {
        hid_t file_id = open_hdf5(filename);
        ssf_manager::write_data(file_id);
        // write_data(file_id);
        H5Fclose(file_id);
    }


};
