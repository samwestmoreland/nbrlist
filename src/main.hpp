#ifndef MAIN_H_
#define MAIN_H_

/* function prototypes for main.cpp */
int set_default_parameters();
int calculate_interactions();
vec_t get_uc_dimensions_from_zr_content(double zrconcentration);
int generate_domain_wall_system();

int generate_large_system(std::vector<pair_t>& uc_interactions,
                          std::vector<atom_t>& unitcell,
                          std::vector < std::vector < std::vector < std::vector <atom_t> > > >& system,
                          vec_t ucd,
                          int n_tracked_cells,
                          int n_materials,
                          int system_dimension);

#endif /* MAIN_H_ */
