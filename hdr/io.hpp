#ifndef IO_H_
#define IO_H_

int output_materials(std::vector<material_t>& materials, std::vector<int>& material_specific_atom_count, int n_atoms);

parameter_t parse_input (std::string const& inputfile);

void array_to_rasmol(std::vector<atom_t> array, std::string const& arrayname);

#endif /* IO_H_ */
