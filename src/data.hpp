#ifndef DATA_H_
#define DATA_H_

#include<vector>

/* arrays (global) */
extern std::vector<atom_t> unitcell;
extern std::vector<atom_t> expanded_uc;
extern std::vector<atom_t> supercell;
extern std::vector<atom_t> expanded_sc;
extern std::vector<std::vector<std::vector<std::vector<atom_t> > > > domainwallsystem;
extern std::vector<int_t> uc_interactions;
extern std::vector<material_t> materials;
extern std::vector<int> material_specific_atom_count;
extern vec_t ucd;
extern vec_t exp_ucd;

/* output file stream */
extern std::ofstream outfile;

#endif /* DATA_H_ */
