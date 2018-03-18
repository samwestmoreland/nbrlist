#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"

parameter_t parse_input (std::string const& inputfile) {

   parameter_t system;

   std::ifstream fin;
   fin.open(inputfile.c_str());

   /* exit if no input file found */
   if (!fin.good()) {
       std::cout << "no input file found. program exiting." << std::endl;
       exit(EXIT_FAILURE);
   }
   else
      std::cout << "reading parameters from file '" << inputfile << "'" << std::endl;

   std::string line;

   /* read each line of file */
   while (std::getline(fin, line))
   {
      /* skip comments or blank lines */
      if ((line[0]=='#')||line.empty()) continue;

      /* remove spaces from line */
      line.erase(remove(line.begin(), line.end(), ' '), line.end());

      /* identify delimiters and initialise key string */
      const char* eq = "=";
      const char* cm = ",";

      std::string key = "";
      std::string val = "";

      /* point in line where key ends */
      int endvar;

      /* loop through characters in line to identify key */
      for (int i=0; i<line.length(); ++i) {

         if (line.at(i) != *eq) key.push_back(line.at(i));
         else {
            endvar = i+1;
            break;
         }
      }

      /* loop through characters after equals sign */
      for (int i=endvar; i<line.length(); ++i) val.push_back(line.at(i));

      if (key == "cutoff") system.rcut = stof(val);

      else if (key == "tt_factor") system.tt_factor = stof(val);
      else if (key == "rt_factor") system.rt_factor = stof(val);

      else if (key == "tracking") {
         if (val == "false") system.tracking = false;
         else if (val == "true") system.tracking = true;
         else std::cout << "input parse error: i don't know what " << "\"" << val << "\" means\n";
      }

      else if (key == "material") {
         system.material = val;
         system.material_int = convert_material_string_to_integer(val);
      }

      else if (key == "config") {
         system.config = val;
      }

      else if (key == "zrconcentration") {
         system.zrconcentration = stof(val);
         if (system.zrconcentration != 0) system.zrdoping = true;
      }

      else if (key == "domainwall") {
         if (val == "true") system.domainwall = true;
         else if (val == "false") system.domainwall = false;
         else {
            std::cout << "not sure what '" << val << "' means. exiting\n";
            exit(EXIT_FAILURE);
         }
      }

      else if (key == "dwsystemdimensionx") {
         system.domainwall = true;
         system.dw_dim.x = stof(val);
      }

      else if (key == "dwsystemdimensiony") {
         system.domainwall = true;
         system.dw_dim.y = stof(val);
      }

      else if (key == "dwsystemdimensionz") {
         system.domainwall = true;
         system.dw_dim.z = stof(val);
      }

      else {
         std::cout << "input parse error: i don't know what " << "\"" << key << "\" means\n";
         exit(EXIT_FAILURE);
      }

   }

   return system;
}

int output_materials(std::vector<material_t>& materials, std::vector<int>& material_specific_atom_count, int n_atoms) {

    for (int i=0; i<materials.size(); ++i)
        std::cout << "\t"
                  << materials[i].name << "\t"
                  << materials[i].id << "\t"
                  << material_specific_atom_count[i] << "\t"
                  << float(material_specific_atom_count[i])/n_atoms*100. << "%\n";

    if (materials[materials.size()-1].name == "Zr") {
       std::cout << std::endl
                 << "zr concentration: "
                 << float(material_specific_atom_count[materials.size()-1])/(float(material_specific_atom_count[0])+float(material_specific_atom_count[materials.size()-1]))
                 << std::endl;
    }

    return EXIT_SUCCESS;
}

void array_to_rasmol(std::vector<atom_t> array, std::string const& arrayname) {

   std::string filename = arrayname + ".xyz";
   std::ofstream rasmol (filename.c_str());

   /* if file open */
   if (rasmol) {

      rasmol << array.size() << "\n\n";

      for (int i=0; i<array.size(); ++i)
         rasmol
            << array[i].element << "\t"
            << array[i].pos.x   << "\t"
            << array[i].pos.y   << "\t"
            << array[i].pos.z   << "\n";

      rasmol.close();
      std::cout << filename << " file generated\n";
   }
}
