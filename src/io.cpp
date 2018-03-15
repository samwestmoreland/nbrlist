#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"

int initialise_material(int material_int, std::string const& material, double zr_content, std::string config) {

   std::cout << "\ninitialising material \'" << material << "\'\n";

   switch(material_int) {

      case 1 :    /* bccfe */
         ucd.x = 2.856;
         ucd.y = 2.856;
         ucd.z = 2.856;
         break;

      case 2 :    /* ndfeb */
         ucd.x = 8.8;
         ucd.y = 8.8;
         ucd.z = 12.2;
         break;

      case 3 :    /* ndfe12 */
         ucd.x = 8.574;
         ucd.y = 8.574;
         ucd.z = 4.907;
         break;

      case 4 :    /* smfe12 */
         ucd.x = 8.497;
         ucd.y = 8.497;
         ucd.z = 4.687;
         break;

      case 5 :    /* smzrfe12 */
         ucd = get_uc_dimensions_from_zr_content(zr_content);
         break;

      case 6 :    /* interface */
         ucd.x = 26.181;
         ucd.y = 26.181;
         ucd.z = 0;
         break;

      case 7 :    /* interface_mirror */
         ucd.x = 26.181;
         ucd.y = 26.181;
         ucd.z = 0;
         break;
   }

   std::cout << std::endl;
   std::string filename = generate_filename(material, config);

   std::ifstream infile (filename.c_str());

   // std::string filename;
   // std::string config_str = std::to_string(config);

   // if (material == "smzrfe12")
   //     filename = "coordinates/smzrfe12/" + material + config_str + ".coords";
   // else filename = "coordinates/" + material + ".coords";

   /* check if file opened */
   if (!(infile.is_open())) {
      std::cerr   << "couldn't find coordinate file \'" << filename << "\'" << std::endl;
      std::cerr   << "exiting" << std::endl;
      std::exit(EXIT_FAILURE);
   }

   std::cout
      << "unit cell dimensions: "
      << ucd.x << ", "
      << ucd.y << ", "
      << ucd.z << "\n";

   atom_t temp;
   int atom_count = 0;

   std::vector<int> material_specific_atom_count;

   std::cout << "reading in unit cell coordinates from " << filename << std::endl;

   while (infile
       >> temp.element
       >> temp.pos.x
       >> temp.pos.y
       >> temp.pos.z) {

            /* assign atom id */
            temp.aid = atom_count;
            atom_count ++;

            /* assign some dummy variables for unneeded struct elements */
            temp.gid = 0;
            temp.hcat = 0;
            temp.uc.x = 0;
            temp.uc.y = 0;
            temp.uc.z = 0;

            /* determine material id */
            material_t temp_mat;
            temp_mat.name = temp.element;
            temp_mat.id = determine_material_id(temp_mat.name, materials);
            temp.mat = temp_mat.id;

            /* make sure material_specific_atom_count has correct number of elements */
            if (temp_mat.id+1 > material_specific_atom_count.size())
               material_specific_atom_count.push_back(0);

            /* add one to the material atom counter */
            material_specific_atom_count[temp_mat.id] ++;

            if (temp.element == "Fe8i" || temp.element == "Fe8j" || temp.element == "Fe8f" || temp.element == "Fe")
               temp.fe = true;
            else temp.fe = false;

            /* scale atom coordinates */
            temp.pos = temp.pos * ucd;

            unitcell.push_back(temp);
      }

   std::cout << "atoms read in: " << unitcell.size() << std::endl;

   std::cout << "number of materials: " << materials.size() << "\n\n\t";

   /* output names of materials */
   output_materials(materials, material_specific_atom_count, unitcell.size());

   /*
      output atom informaton to ucf
   */

   outfile.open("output.ucf");

   outfile
      << "# Unit cell size:\n"
      << ucd.x << "\t"
      << ucd.y << "\t"
      << ucd.z << "\n"
      << "# Unit cell vectors\n"
      << "1.0 0.0 0.0\n"
      << "0.0 1.0 0.0\n"
      << "0.0 0.0 1.0\n"
      << "# Atoms num, id cx cy cz mat lc hc\n"
      << unitcell.size() << std::endl;

   for (int i=0; i<unitcell.size(); ++i) {

      outfile
         << unitcell[i].aid << "\t"
         << unitcell[i].pos.x/ucd.x << "\t"
         << unitcell[i].pos.y/ucd.y << "\t"
         << unitcell[i].pos.z/ucd.z << "\t"
         << unitcell[i].mat << "\t"
         << "0\t0\n";
   }

   return EXIT_SUCCESS;
}

parameter_t parse_input (std::string const& inputfile)
{
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
        std::cout << materials[i].name << "\t" << materials[i].id << "\t" << material_specific_atom_count[i] << "\t" << float(material_specific_atom_count[i])/n_atoms*100. << "%\n\t";

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
      std::cout << filename << " file generated\n\n";
   }
}
