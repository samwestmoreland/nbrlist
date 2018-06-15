#ifndef INITIALISE_H_
#define INITIALISE_H_

/* declarations of functions defined in initialise.cpp */
int determine_element_id(std::string const& element);
int get_unitcell_dimensions();
int read_coordinatefile();
std::string generate_filename(std::string const& material_string);

#endif /* INITIALISE_H_ */
