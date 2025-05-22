#pragma once 
#include <cstddef>
#include <string>
#include <vector>
#include <sstream>
#include <lattice_IO.hpp>
#include <argparse.hpp>


template<typename T>
inline std::string comma_separate(const char* pname, std::vector<T> v){
    std::ostringstream oss;
    oss<<pname;
    char delim='=';
    for (auto n : v){
        oss << delim << n;
        delim=','; 
    }
    oss << ";";
    return oss.str();
}


inline void parse_supercell_spec(imat33_t& supercell_spec, 
        const argparse::ArgumentParser& prog){


    int L = prog.get<size_t>("L");

    std::vector<int> Z1 = {L, 0, 0};
    std::vector<int> Z2 = {0, L, 0};
    std::vector<int> Z3 = {0, 0, L};


    for (int row=0; row<3; row++){
        supercell_spec(row,0) = Z1[row];
        supercell_spec(row,1) = Z2[row];
        supercell_spec(row,2) = Z3[row];
    }
}



inline uint64_t get_seed(const argparse::ArgumentParser& prog){ 
    std::string seed_s = prog.get<std::string>("--seed");
    std::stringstream ss;
    ss << std::hex << seed_s;
    uint64_t seed;
    ss >> seed;
    return seed;
}

