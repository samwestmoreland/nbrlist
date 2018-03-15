#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"
#include "./initialise.hpp"

int initialise_material(int material_int, std::string const& material, double zrconcentration) {

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
         ucd = get_uc_dimensions_from_zr_content(zrconcentration);
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

   std::string filename = generate_filename(material, zrconcentration);
   std::ifstream infile (filename.c_str());

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

std::string generate_filename(std::string const& material_string, double zrconcentration) {

   std::string filename;

   /* convert double zrconcentration to string */
   std::ostringstream zr_strs;
   zr_strs << zrconcentration;
   std::string zr_string = zr_strs.str();
   /**/

   if (material_string == "smzrfe12")
      filename = "./coordinates/smzrfe12/config" + zr_string + ".coords";
   else
      filename = "./coordinates/" + material_string + ".coords";

   return filename;
}
