#include <fstream>

#include "./classes.hpp"
#include "./data.hpp"

/* global atom arrays */
std::vector<atom_t> unitcell;
std::vector<atom_t> supercell;
std::vector<std::vector<std::vector<std::vector<atom_t> > > > domainwallsystem;

/* unit cell interaction array */
std::vector<pair_t> uc_interactions;

/* element arrays */
std::vector<std::string> elements;
std::vector<int> element_specific_atom_count;

/* global parameter array */
parameter_t sys;

/* material paramters */
material_t mat;

/* global output file stream */
std::ofstream outfile;
