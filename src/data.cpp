#include <fstream>

#include "../hdr/classes.hpp"
#include "../hdr/data.hpp"

/* arrays (global) */
std::vector<atom_t> unitcell;
std::vector<atom_t> supercell;
std::vector<int_t> uc_interactions;
std::vector<material_t> materials;
vec_t ucd;

/* output file stream */
std::ofstream outfile;

