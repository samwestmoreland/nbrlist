#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"

void parse_input (std::string const& inputfile) {

   std::ifstream fin;
   fin.open(inputfile.c_str());

   /* exit if no input file found */
   if (!fin.good()) {
       std::cout << "no input file found. program exiting." << std::endl;
       exit(EXIT_FAILURE);
   }

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

      if (key == "tmtmcutoff") sys.tmtmrcut = stof(val);
      else if (key == "retmcutoff") sys.retmrcut = stof(val);

      else if (key == "tt_factor") sys.tt_factor = stof(val);
      else if (key == "rt_factor") sys.rt_factor = stof(val);

      else if (key == "retmexchangeconstant") sys.rt_exchange_constant = stof(val);

      else if (key == "tracking") {
         if (val == "false") sys.tracking = false;
         else if (val == "true") sys.tracking = true;
         else std::cout << "input parse error: i don't know what " << "\"" << val << "\" means\n";
      }

      else if (key == "material") {
         sys.material = val;
         sys.material_int = convert_material_string_to_integer(val);
      }

      else if (key == "zrconcentration") {
         sys.zrconcentration = stof(val);
         if (sys.zrconcentration != 0) sys.zrdoping = true;
      }

      else if (key == "domainwall") {
         if (val == "true") sys.domainwall = true;
         else if (val == "false") sys.domainwall = false;
         else {
            std::cout << "not sure what '" << val << "' means. exiting\n";
            exit(EXIT_FAILURE);
         }
      }

      else if (key == "centrepin") {

         if (val == "true") sys.centrepin = true;
         else if (val == "false") sys.centrepin = false;
         else {
            std::cout << "not sure what '" << val << "' means. exiting\n";
            exit(EXIT_FAILURE);
         }
      }

      else if (key == "dwsystemdimensionx") {
         sys.dw_dim.x = stof(val);
      }

      else if (key == "dwsystemdimensiony") {
         sys.dw_dim.y = stof(val);
      }

      else if (key == "dwsystemdimensionz") {
         sys.dw_dim.z = stof(val);
      }

      else if (key == "ticoncentration") {
         sys.ticoncentration = stof(val);
      }

      else {
         std::cout << "input file parse error: i don't know what " << "\"" << key << "\" means\n";
         exit(EXIT_FAILURE);
      }

   }

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
