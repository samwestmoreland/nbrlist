#ifndef DATA_H_
#define DATA_H_

#include <vector>

#include "./classes.hpp"

/* global atom arrays */
extern std::vector<atom_t> unitcell;
extern std::vector<atom_t> supercell;
extern std::vector<std::vector<std::vector<std::vector<atom_t> > > > domainwallsystem;

/* unit cell interaction array */
extern std::vector<pair_t> uc_interactions;

/* element arrays */
extern std::vector<std::string> elements;
extern std::vector<int> element_specific_atom_count;

extern vec_t ucd;
extern parameter_t sys;
extern material_t mat;

/* output file stream */
extern std::ofstream outfile;

#endif /* DATA_H_ */
