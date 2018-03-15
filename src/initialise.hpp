#ifndef INIT_H_
#define INIT_H_

/* prototypes for functions defined in initialise.cpp */
int initialise_material(int material_int, std::string const& material, double zrconcentration);
std::string generate_filename(std::string const& material_string, double zrconcentration);

#endif /* INIT_H_ */
