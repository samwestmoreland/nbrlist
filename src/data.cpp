#include <fstream>

#include "./classes.hpp"
#include "./data.hpp"

/* arrays (global) */
std::vector<atom_t> unitcell;
std::vector<atom_t> supercell;
std::vector<std::vector<std::vector<std::vector<atom_t> > > > domainwallsystem;
std::vector<int_t> uc_interactions;
std::vector<material_t> materials;
std::vector<int> material_specific_atom_count;
vec_t ucd;
parameter_t sys;

/* global output file stream */
std::ofstream outfile;
