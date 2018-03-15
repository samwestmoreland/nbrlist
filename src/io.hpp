#ifndef IO_H_
#define IO_H_

int convert_material_string_to_integer(std::string const& material);
vec_t get_uc_dimensions_from_zr_content(double zr_content);
int determine_material_id(std::string const& in_material, std::vector<material_t>& materials);
parameter_t parse_input (std::string const& inputfile);
void array_to_rasmol(std::vector<atom_t> array, std::string const& arrayname);
int output_materials(std::vector<material_t>& materials, std::vector<int>& material_specific_atom_count, int n_atoms);

#endif /* IO_H_ */
