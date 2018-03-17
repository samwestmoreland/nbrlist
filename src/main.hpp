#ifndef MAIN_H_
#define MAIN_H_

/* function prototypes for main.cpp */
int convert_material_string_to_integer(std::string const& material);
double calculate_rij(vec_t& i, vec_t& j);     /* distance calculation */
double jij_ndfeb(std::string const& i_type, std::string const& j_type, double rij);
int calculate_interactions(int exchange_fn, double tt_factor, double rt_factor, double rcut, double zrconcentration);
vec_t get_uc_dimensions_from_zr_content(double zrconcentration);
int expand_unitcell_and_substitute_zr_atoms(double zrconcentration);
int generate_domain_wall_system();

int generate_large_system(std::vector<int_t>& uc_interactions,
                          std::vector<atom_t>& unitcell,
                          std::vector < std::vector < std::vector < std::vector <atom_t> > > >& system,
                          vec_t ucd,
                          int n_tracked_cells,
                          int n_materials,
                          int system_dimension);

double calculate_jij(std::string const& i_type,
                     std::string const& j_type,
                     double rij,
                     double tt_factor,
                     double rt_factor,
                     int exchange_fn);

#endif /* MAIN_H_ */
