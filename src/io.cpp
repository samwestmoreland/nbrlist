#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"

int output_header() {

   std::cout << std::endl;
   std::cout << "***************************************************\n\n";
   std::cout << "                   Nbrlist\n\n";
   std::cout << "       Compiled on " << __DATE__ << " at " << __TIME__ << "\n\n";
   std::cout << "***************************************************\n";

   return EXIT_SUCCESS;
}

int read_inputfile() {

   std::cout << "\nreading input file...";
   std::ifstream inputfile ("nbrlist.in");

   /* exit if no input file found */
   if (!inputfile.good()) {
       std::cout << "no input file found. program exiting." << std::endl;
       exit(EXIT_FAILURE);
   }

   std::string line;

   /* read each line of file */
   while (std::getline(inputfile, line))
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

      if (key == "nearestneighbour") sys.nn = true;

      else if (key == "tmtmcutoff") sys.rcut_tt = stof(val);
      else if (key == "retmcutoff") sys.rcut_rt = stof(val);

      else if (key == "ttfactor") sys.tt_factor = stof(val);
      else if (key == "rtfactor") sys.rt_factor = stof(val);

      else if (key == "rtshell") {
         sys.rt_shell = stoi(val);
         sys.rt_shell_cut = true;
      }

      else if (key == "rtcutoff") {
         sys.rcut_rt = stof(val);
         sys.rt_shell_cut = false;
      }

      else if (key == "retmexchangeconstant") sys.rt_constant = stof(val);
      else if (key == "tmtmexchangeconstant") sys.tt_constant = stof(val);

      else if (key == "tracking") {
         if (val == "false") sys.tracking = false;
         else if (val == "true") sys.tracking = true;
         else std::cout << "input parse error: i don't know what " << "\"" << val << "\" means\n";
      }

      else if (key == "material") {
         mat.name = val;
         if (mat.id() == 0) {
            std::cout << "material \'" << val << "\' not recognised. exiting.\n";
            exit(EXIT_FAILURE);
         }
      }

      else if (key == "zrconcentration") {
         sys.zr_concentration = stof(val);
      }

      else if (key == "ticoncentration") {
         sys.ti = stof(val);
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

      else {
         std::cout << "input file parse error: i don't know what " << "\"" << key << "\" means\n";
         exit(EXIT_FAILURE);
      }

   }

   /* output read parameters */
   std::cout << "\nparameters read from input file:\n\n";
   std::cout << "\tmaterial: " << mat.name << std::endl;
   std::cout << "\ttm-tm cut-off radius: " << sys.rcut_tt << " A" << std::endl;

   if (sys.rt_shell_cut) std::cout << "\tre-tm cut-off shell: " << sys.rt_shell << std::endl;
   else std::cout << "\tre-tm cut-off radius: " << sys.rcut_rt << " A" << std::endl;

   std::cout << std::endl;
   std::cout << "\tzr concentration: " << sys.zr_concentration << std::endl;
   std::cout << "\tti concentration: " << sys.ti << std::endl;
   std::cout << std::endl;
   std::cout << "\tt-t exchange factor = " << sys.tt_factor << std::endl;
   std::cout << "\tr-t exchange constant = " << sys.rt_constant << std::endl;

   return EXIT_SUCCESS;

}

int output_elements() {

    for (int i=0; i<elements.size(); ++i)
        std::cout
           << "\t" << elements[i]
           << "\t" << element_specific_atom_count[i]
           << "\t" << float(element_specific_atom_count[i])/unitcell.size()*100. << "%\n";

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
