#ifndef INIT_H_
#define INIT_H_

int initialise_material(int material_int, std::string const& material, double zrconcentration);
std::string generate_filename(std::string const& material_string, double zrconcentration);

#endif /* INIT_H_ */
