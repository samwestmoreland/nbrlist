#ifndef MAIN_H_
#define MAIN_H_

/* function prototypes */
int convert_material_string_to_integer(std::string const& material);

double calculate_rij(vec_t& i, vec_t& j);     /* distance calculation */

double calculate_jij(std::string const& i_type,
                     std::string const& j_type,
                     double rij,
                     double tt_factor,
                     double rt_factor,
                     int exchange_fn);

double jij_ndfeb(std::string const& i_type,
                 std::string const& j_type,
                 double rij);

int initialise_material(int material_int, std::string const& material, double zr_content, std::string config);

int calculate_interactions(int exchange_fn,
                           double tt_factor,
                           double rt_factor,
                           double rcut);



int generate_large_system(std::vector<int_t>& uc_interactions,
                          std::vector<atom_t>& unitcell,
                          std::vector < std::vector < std::vector < std::vector <atom_t> > > >& system, vec_t ucd,
                          int n_tracked_cells,
                          int n_materials,
                          int system_dimension);

#endif /* MAIN_H_ */
