#ifndef INITIALISE_H_
#define INITIALISE_H_

/* declarations of functions defined in initialise.cpp */
int determine_element_id(std::string const& element);
int determine_element_id_for_interface_system(std::string const& element, double zpos);
int get_unitcell_dimensions();
int read_coordinatefile();
std::string generate_filename(std::string const& material_string);
int neighbour_routines();
int identify_re_neighbours();
int identify_tm_neighbours();

#endif /* INITIALISE_H_ */
