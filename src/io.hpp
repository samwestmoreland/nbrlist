#ifndef IO_H_
#define IO_H_

int read_inputfile();
int output_header();
void array_to_rasmol(std::vector<atom_t> array, std::string const& arrayname);
int output_elements();

#endif /* IO_H_ */
