#include "MC_routines.hpp"
#include "format_bits.hpp"
#include "observables.hpp"
#include <cstdio>
#include <ostream>
#include <vector>


using namespace std;


void declare_args(argparse::ArgumentParser& prog){

    prog.add_argument("L")
        .help("Linear system size")
        .scan<'u', size_t>();

    prog.add_argument("n_sweeps")
        .help("Number of MC sweeps")
        .scan<'u', size_t>();

    prog.add_argument("n_loop")
        .help("Number of MC sweeps")
        .scan<'u', size_t>();

    prog.add_argument("n_local")
        .help("Number of MC sweeps")
        .scan<'u', size_t>();

    prog.add_argument("--Jzz")
        .help("The six couplings (between sls 01, 02, 03, 23, 13, 12) is the coupling sl1<->sl2)")
        .nargs(6)
        .required()
        .scan<'g', double>();

    prog.add_argument("--output_dir", "-o")
        .help("Path to output")
        .required();

    prog.add_argument("--temperature", "-T")
        .help("Temperature (can be 0)")
        .scan<'g', double>()
        .required();

    prog.add_argument("--sweep_strategy", "-s")
        .help("Sweep strategy: (s)equential or (r)andomised")
        .choices("s", "r")
        .default_value("s");

    prog.add_argument("--max_loop_len", "-l")
        .help("Largest loop allowed in the loop MC algorithm")
        .scan<'u', size_t>()
        .default_value(std::numeric_limits<size_t>::max());

    prog.add_argument("--seed")
        .help("64-bit int to seed the RNG")
        .default_value("0");

}



inline void hash_parameters(std::stringstream &name,
                            const argparse::ArgumentParser &prog) {

    name << "L=" << prog.get<size_t>("L") << "&";
    name << "n_sweeps=" << prog.get<size_t>("n_sweeps") << "&";
    name << "n_loop=" << prog.get<size_t>("n_loop") << "&";
    name << "n_local=" << prog.get<size_t>("n_local") << "&";

    name << "T=" << prog.get<double>("--temperature") << "&";
    name << "seed=" << prog.get("--seed") << "&";
}




inline MC_settings parse_MC_settings(argparse::ArgumentParser& prog){
    MC_settings settings;
    settings.T = prog.get<double>("--temperature");
    auto sweep_strat = prog.get("--sweep_strategy")[0];
    if (sweep_strat == 's'){
        settings.sweep_strategy  = MC_settings::sweep::sequential;
    } else if (sweep_strat == 'r'){
        settings.sweep_strategy  = MC_settings::sweep::random;
    }
    settings.max_loop_len = prog.get<size_t>("--max_loop_len");

    settings.n_sweeps = prog.get<size_t>("n_sweeps");
    settings.sweep_n_local = prog.get<size_t>("n_local");
    settings.sweep_n_loop = prog.get<size_t>("n_loop");
    return settings;
}

////////////////////////////////////////////////////////////////////////////////
///////   MAIN PROGRAM                   ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int main (int argc, const char *argv[]) {
    // Interface
    argparse::ArgumentParser prog(argv[0]);
    
    declare_args(prog);

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err){
        cerr << err.what() << endl;
        cerr << prog;
        std::exit(1);
    }
    //////////////////////////////////////////////////////// 
    /// End program argument definitions
    ///
    std::string output_dir = prog.get("output_dir");
    std::filesystem::path outpath(output_dir);
    if (! filesystem::exists(outpath) ){
        throw std::runtime_error("Cannot open output_dir");
    }


    std::stringstream name; // accumulates hashed options
    hash_parameters(name, prog);

    // parse L1 L2 L3
    imat33_t supercell_spec;
    parse_supercell_spec(supercell_spec, prog);
    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    // option setting done: build the lattice
    auto ssf_path = outpath/(name.str()+".ssf.h5");

    const auto spec = PrimitiveSpecifiers::DiamondSpec();

    auto settings = parse_MC_settings(prog);


    /*
     * Here, we intitalize all containers needed.
     *
     */
    Lattice lat(spec, supercell_spec);
    
    auto Jzz = prog.get<std::vector<double>>("--Jzz");
    // 01 02 03 23 31 12
    static const int Jzz_index[4][4] = {
        {-1,  0,  1,  2},
        { 0, -1,  5,  4},
        { 1,  5, -1,  3},
        { 2,  4,  3, -1}
    };

    setup_couplings(lat, [Jzz](size_t i,size_t j){return Jzz[Jzz_index[i][j]];});

    rng_t gen(get_seed(prog));

    cout << "Seed " <<get_seed(prog) <<std::endl;

    initialise_spins(lat, gen);


    std::vector<Spin*> spins;
    for (auto [_, s] : lat.links){
        spins.push_back(s);
    }
    
    auto ssf_manager = SSF_manager(lat);

    FILE* dbf = fopen("debug.txt", "w");

    for (size_t n=0; n<settings.n_sweeps; n++){

        for (auto& [i, t] : lat.points ){
            fprintf(dbf, "%2d [", i);
            for (auto& [l, _] : t->coboundary) {
                auto s = static_cast<Spin*>(l);
                fprintf(dbf, "%d",(s->sz+1)/2);
            }
            fprintf(dbf, "] ");
        }
        fprintf(dbf,"\n");


        size_t success_local=0;
        size_t success_loop=0;
        for (size_t i=0; i<settings.sweep_n_local; i++){
            success_local += mc_sweep_local(spins, gen, settings);
        }
        for (size_t i=0; i<settings.sweep_n_loop; i++){
            success_loop += mc_sweep_loop(spins, gen, settings);
        }
        ssf_manager.store();
        
        
        printf("Acceptance rate: Local %.6f Loop %.6f\r",
                success_local*1.0/4/lat.num_primitive/settings.sweep_n_local,
                success_loop*1.0/4/lat.num_primitive/settings.sweep_n_loop 
              );
        cout<<std::flush;
    }
    cout <<std::endl;
    fclose(dbf);

    ssf_manager.save_to_hdf5(ssf_path);

    return 0;
}

