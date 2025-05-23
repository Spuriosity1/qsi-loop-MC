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

/*
template<typename Derived> 
struct Observable {
    Observable(){};

    void store(const Lattice& lat) {
        static_cast<Derived*>(this)->store(lat);
    }

    void write_data(hid_t file_id) {
        static_cast<Derived*>(this)->write_data(file_id);
    }

};
*/

template <typename T>
concept ObservableType = requires(const Lattice &lat, T t, hid_t hid) {
    T{lat}; // requires a constructor T(const Lattice&)
    { t.store(lat) } -> std::same_as<void>;
    { t.write_data(hid) } -> std::same_as<void>;
};

hid_t open_hdf5(const std::string &filename);
void write_lat_data( hid_t file_id, const Lattice& lat);
void write_n_samples(hid_t file_id, size_t n_samples);


template<ObservableType... Obs>
class StatsManager : public Obs... {
public:
    StatsManager(const Lattice& _lat) : Obs(_lat)..., lat(_lat) {}

    void store() {
        n_samples++;
        (Obs::store(lat), ...);  // fold expression
    }

    void save(const std::string& filename) {
        hid_t file_id = open_hdf5(filename); 
        write_lat_data(file_id, lat);
        write_n_samples(file_id, n_samples);
        
        (Obs::write_data(file_id), ...);
        H5Fclose(file_id);
    }
private:
    const Lattice& lat;
    uint64_t n_samples = 0;
};


// class responsible for calculating \sum_i \sigma_i e^{iq \cdot r}
struct SSFObservable {
    SSFObservable(const Lattice& lat); 
    ~SSFObservable();
    void store(const Lattice& lat);
    void write_data(hid_t file_id);

    
private:

    fftw_plan plans[4];
    // temporary storage
    fftw_complex *fft_out;
    double *spin_field[4];
    // the accumulator for <Sz(q)Sz(-q)> (note it is real)
    double *S_q_acc; // indexed by sublattice
    hsize_t S_q_acc_dims[3];

};



// adds on raw <Sz> measurements
struct SzObservable {
    SzObservable(const Lattice& lat) {
        Sz_acc = new double[lat.num_primitive * 4];
        std::memset(Sz_acc, 0, sizeof(double) * lat.num_primitive * 4);

        const auto L = lat.size();
        Sz_dims[0]= 4;
        for (int d=0; d<3; d++){
            Sz_dims[d+1] = static_cast<hsize_t>(L[2]);
        }
    }
    ~SzObservable(){ delete[] Sz_acc; }

    void store(const Lattice& lat){
        for (auto& [J, l] : lat.links){
            auto s = static_cast<Spin*>(l);
            Sz_acc[J] += s->sz;
        }
    }
    void write_data(hid_t file_id);


private:
    double *Sz_acc;

    hsize_t Sz_dims[4]; 

};

