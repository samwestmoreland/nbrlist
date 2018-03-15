#ifndef DATA_H_
#define DATA_H_

#include<vector>

/* arrays (global) */
extern std::vector<atom_t> unitcell;
extern std::vector<atom_t> supercell;
extern std::vector<int_t> uc_interactions;
extern std::vector<material_t> materials;
extern vec_t ucd;

/* output file stream */
extern std::ofstream outfile;

#endif /* DATA_H_ */
