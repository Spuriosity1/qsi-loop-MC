#include "observables.hpp"
#include "struct.hpp"

/// Base class, geometry storage

hid_t open_hdf5(const std::string &filename) {
    H5Eset_auto(H5E_DEFAULT, (H5E_auto2_t) H5Eprint, stderr);
    return H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

void write_n_samples(hid_t file_id, size_t n_samples){

    // Save number of samples
    {
        hsize_t dims[1] = {1};
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t dset = H5Dcreate2(file_id, "n_samples", H5T_NATIVE_UINT64, space,
                                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          &n_samples);
        H5Dclose(dset);
        H5Sclose(space);
    }

}

void write_lat_data( hid_t file_id, const Lattice& lat) {

    // Save lattice size
    {
      hsize_t dims[1] = {3};
      long size[3] = {lat.size()[0], lat.size()[1], lat.size()[2]};
      hid_t space = H5Screate_simple(1, dims, nullptr);
      hid_t dset = H5Dcreate2(file_id, "lattice_size", H5T_NATIVE_INT, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, size);
      H5Dclose(dset);
      H5Sclose(space);
    }

    // Save basis vectors
    {
        hsize_t dims[2] = {3, 3};
        double basis[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                basis[i][j] = convert_rational<double>(
                        lat.cell_vectors(i, j)); // Eigen column-major

        hid_t space = H5Screate_simple(2, dims, nullptr);
        hid_t dset = H5Dcreate2(file_id, "basis_vectors", H5T_NATIVE_DOUBLE, space,
                                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, basis);
        H5Dclose(dset);
        H5Sclose(space);
    }


    // Save primitive basis vectors
    {
        hsize_t dims[2] = {3, 3};
        double basis[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                basis[i][j] = convert_rational<double>(
                        lat.primitive_spec.lattice_vectors(i, j)); // Eigen column-major

        hid_t space = H5Screate_simple(2, dims, nullptr);
        hid_t dset = H5Dcreate2(file_id, "primitive_basis_vectors", H5T_NATIVE_DOUBLE, space,
                                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, basis);
        H5Dclose(dset);
        H5Sclose(space);
    }

    // Save sublattice offsets
    {

        hsize_t n_sl = lat.primitive_spec.num_link_sl();
        hsize_t dims[2] = {n_sl, 3};
        std::vector<double> flat(n_sl * 3);
        for (size_t sl = 0; sl < n_sl; ++sl)
            for (int j = 0; j < 3; ++j)
                flat[sl * 3 + j] = convert_rational<double>(
                        lat.primitive_spec.link_no(sl).position[j]);

        hid_t space = H5Screate_simple(2, dims, nullptr);
        hid_t dset = H5Dcreate2(file_id, "sublattice_offsets", H5T_NATIVE_DOUBLE,
                                                        space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          flat.data());
        H5Dclose(dset);
        H5Sclose(space);
    }
}


// SSF storage

SSFObservable::SSFObservable(const Lattice &lat) {
    auto Lvec = lat.size();
    assert(Lvec[0] == Lvec[1] && Lvec[1] == Lvec[2]); // otherwise it's too hard
	auto L = Lvec[0];	
	
    fft_out = fftw_alloc_complex(L * L * (L / 2 + 1));

	S_q_acc_dims[0] = 4; // num sublattices
	S_q_acc_dims[1] = L+1;
	S_q_acc_dims[2] = L/2+1;

    for (int sl = 0; sl < 4; sl++) {
        spin_field[sl] = fftw_alloc_real(lat.num_primitive);
        plans[sl] = fftw_plan_dft_r2c_3d(L, L, L, spin_field[sl], fft_out,
                                         FFTW_ESTIMATE);
    }

    auto N = S_q_acc_dims[0] * S_q_acc_dims[1] * S_q_acc_dims[2];
    S_q_acc = new double[ N ];
    std::memset(S_q_acc, 0, sizeof(double) * N);
}






/*
ssf_manager::ssf_manager(const Lattice &_lat) : 
    observable_manager(_lat) {
    auto Lvec = lat.size();
    assert(Lvec[0] == Lvec[1] && Lvec[1] == Lvec[2]); // otherwise it's too hard
	auto L = Lvec[0];	
	
    fft_out = fftw_alloc_complex(L * L * (L / 2 + 1));

	S_q_acc_dims[0] = 4; // num sublattices
	S_q_acc_dims[1] = L+1;
	S_q_acc_dims[2] = L/2+1;

    for (int sl = 0; sl < 4; sl++) {
        spin_field[sl] = fftw_alloc_real(lat.num_primitive);
        plans[sl] = fftw_plan_dft_r2c_3d(L, L, L, spin_field[sl], fft_out,
                                         FFTW_ESTIMATE);
    }

    auto N = S_q_acc_dims[0] * S_q_acc_dims[1] * S_q_acc_dims[2];
    S_q_acc = new double[ N ];
    for (hsize_t i=0; i< N; i++){
        S_q_acc[i] = 0.;
    }
}
*/

SSFObservable::~SSFObservable() {
    for (int sl = 0; sl < 4; sl++) {
        fftw_destroy_plan(plans[sl]);
        fftw_free(spin_field[sl]);
    }
    fftw_free(fft_out);

    delete[] S_q_acc;
}


void SSFObservable::store(const Lattice& lat) { // saves current state of the lattice to	
    auto L = lat.size()[0];
    for (auto [idx, l] : lat.links) {
		// exploits the structure of lat.links' indexing system
		// TODO this should really be provided by LIL itself, since it uses forbidden internal knowledge
        auto sl = idx / lat.num_primitive;
		assert( lat.primitive_spec.sl_of_link(l->position) == sl);
		auto tmp = idx % lat.num_primitive;
		int J[3];
		J[0] = tmp / (L*L);
		J[1] = (tmp - L*L*J[0]) / L;
		J[2] = tmp - L*(L*J[0]+J[1]);
		assert(lat.primitive_spec.lattice_vectors * ipos_t({J[2],J[1],J[0]}) + lat.primitive_spec.link_no(sl).position == l->position);


        assert(sl < 4);
        spin_field[sl][idx % lat.num_primitive] = l->sz;
    }


	

    for (hsize_t sl = 0; sl < S_q_acc_dims[0]; sl++) {
        fftw_execute(plans[sl]);
        // result now available in fft_out
		for (int h=-L/2; h <= L/2; h++){
			for (int l=-L/2; l < L/2; l+=2){ // must jump by 2 to ensure h+l remains even

				auto _i = ((h+l)/2 +L) % L;
				auto _j = ((h+l)/2 +L) % L;
				auto _k = (h+L) % L;

				_k = (_k < L/2 +1) ? _k : L - _k;
				
				// to be correct we need to wrap these
				
				int idx = (_i*L + _j) * (L/2 + 1) + _k;

                double re = fft_out[idx][0];
                double im = fft_out[idx][1];

                double Sq = (re * re + im * im) / (L * L * L);
				
				S_q_acc[sl * S_q_acc_dims[2]*S_q_acc_dims[1] 
                    + (h +L/2)* S_q_acc_dims[2] + (l +L/2)/2 ] += Sq;
			}
		}

    }
}


void SSFObservable::write_data(hid_t file_id) {
    // Save S(q) for each sublattice
    {
        std::string name = "S_q_sublattice";
        //hsize_t dims[2] = {Lx, Nk};
        hid_t space = H5Screate_simple(3, S_q_acc_dims, nullptr);
        hid_t dset = H5Dcreate2(file_id, name.c_str(), H5T_NATIVE_DOUBLE, space,
                                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          S_q_acc);
        H5Dclose(dset);
        H5Sclose(space);
    }

} 

// Sz storage

void SzObservable::write_data(hid_t file_id) {

    // Save S(q) for each sublattice
    {
        std::string name = "Sz_data";
        hid_t space = H5Screate_simple(4, Sz_dims, nullptr);
        hid_t dset = H5Dcreate2(file_id, name.c_str(), H5T_NATIVE_DOUBLE, space,
                                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          Sz_acc);
        H5Dclose(dset);
        H5Sclose(space);
    }

} 
