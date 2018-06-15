#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"
#include "./initialise.hpp"
#include "./main.hpp"

int get_unitcell_dimensions() {

   switch(mat.id()) {

      case 1 :    /* bccfe */
         mat.ucd.x = 2.856;
         mat.ucd.y = 2.856;
         mat.ucd.z = 2.856;
         break;

      case 2 :    /* ndfeb */
         mat.ucd.x = 8.8;
         mat.ucd.y = 8.8;
         mat.ucd.z = 12.2;
         break;

      case 3 :    /* ndfe12 */
//         ucd.x = 8.574;       /*          */
//         ucd.y = 8.574;       /* old code */
//         ucd.z = 4.907;       /*          */

         mat.ucd.x = 8.512;
         mat.ucd.y = 8.512;
         mat.ucd.z = 4.842;     /*  from connor  */
         break;

      case 4 :    /* smfe12 */
         mat.ucd.x = 8.497;
         mat.ucd.y = 8.497;
         mat.ucd.z = 4.687;
         break;

      case 5 :    /* smco12 */
         mat.ucd.x = 8.443068;
         mat.ucd.y = 8.443068;
         mat.ucd.z = 4.799171;    /* from connor */
         break;

      case 6 :    /* interface */
         mat.ucd.x = 26.181;
         mat.ucd.y = 26.181;
         mat.ucd.z = 0;
         break;

   }

   std::cout << std::endl;

   std::cout << "unit cell dimensions:\n\n";
   std::cout << "\tx = " << mat.ucd.x << std::endl;
   std::cout << "\ty = " << mat.ucd.y << std::endl;
   std::cout << "\tz = " << mat.ucd.z << std::endl;

   std::cout << std::endl;

   return EXIT_SUCCESS;

}

int read_coordinatefile() {

   std::string filename = "./coordinates/" + mat.name + ".coords";

   std::ifstream fin (filename.c_str());

   /* check if file opened */
   if (!fin.good()) {
      std::cout << "couldn't open coordinate file \'" << filename << "\'. exiting." << std::endl;
      std::exit(EXIT_FAILURE);
   }

   atom_t tmp;
   int atom_count = 0;

   std::cout << "reading unit cell coordinates from \'" << filename << "\'...";

   while (fin >> tmp.element >> tmp.pos.x >> tmp.pos.y >> tmp.pos.z) {

      /* assign atom id */
      tmp.aid = atom_count;
      atom_count ++;

      /* assign some dummy variables for unneeded struct elements */
      tmp.gid = 0;
      tmp.hcat = 0;
      tmp.uc.x = 0;
      tmp.uc.y = 0;
      tmp.uc.z = 0;

      /* determine element id */
      tmp.mat = determine_element_id(tmp.element);

      /* scale atom coordinates */
      tmp.pos = tmp.pos * mat.ucd;

      unitcell.push_back(tmp);
   }

   std::cout << "done.\n";

   std::cout << "atoms read in: " << unitcell.size() << std::endl;

   std::cout << "number of different elements: " << elements.size() << "\n\n";

   /* output names of elements read */
   output_elements();

   /*
      output atom informaton to ucf
   */

   outfile.open("output.ucf");

   outfile
      << "# Unit cell size:\n"
      << mat.ucd.x << "\t"
      << mat.ucd.y << "\t"
      << mat.ucd.z << "\n"
      << "# Unit cell vectors\n"
      << "1.0 0.0 0.0\n"
      << "0.0 1.0 0.0\n"
      << "0.0 0.0 1.0\n"
      << "# Atoms num, id cx cy cz mat lc hc\n"
      << unitcell.size() << std::endl;

   for (int i=0; i<unitcell.size(); ++i) {

      outfile
         << unitcell[i].aid << "\t"
         << unitcell[i].pos.x/mat.ucd.x << "\t"
         << unitcell[i].pos.y/mat.ucd.y << "\t"
         << unitcell[i].pos.z/mat.ucd.z << "\t"
         << unitcell[i].mat << "\t"
         << "0\t0\n";
   }

   return EXIT_SUCCESS;
}

/* function definition */
std::string generate_filename(std::string const& material_string) {

   std::string filename;

   /* convert double zrconcentration to string */
   std::ostringstream zr_strs;
   zr_strs << sys.zr_concentration;
   std::string zr_string = zr_strs.str();
   /*****/

   filename = "./coordinates/" + material_string + "/" + material_string + ".coords";

   return filename;
}
