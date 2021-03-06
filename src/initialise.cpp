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
         mat.ucd.x = 26.403;
         mat.ucd.y = 26.403;
         mat.ucd.z = 255.93180;
         break;

      case 7 :    /* ndfeti12 */
         if (sys.ti == 0) {
            mat.ucd.x = 17.025013;
            mat.ucd.y = 17.025013;
            mat.ucd.z = 4.842164;
         }

         else if (sys.ti == 1) {
            mat.ucd.x = 17.040534;
            mat.ucd.y = 17.046127;
            mat.ucd.z = 4.8449;
         }

         else if (sys.ti == 2) {
            mat.ucd.x = 17.061906;
            mat.ucd.y = 17.061906;
            mat.ucd.z = 4.847393;
         }

         else if (sys.ti == 3) {
            mat.ucd.x = 17.07819;
            mat.ucd.y = 17.082273;
            mat.ucd.z = 4.849857;
         }

         else if (sys.ti == 4) {
            mat.ucd.x = 17.098785;
            mat.ucd.y = 17.098785;
            mat.ucd.z = 4.85208;
         }

         else {
            std::cout << "no titanium concentration specified. exiting.\n";
            exit(EXIT_FAILURE);
         }

         break;

      case 8 : /* generic */

         mat.ucd.x = 3;
         mat.ucd.y = 3;
         mat.ucd.z = 3;

   }

   std::cout << std::endl;

   std::cout << "unit cell dimensions:\n\n";
   std::cout << "\tx = " << mat.ucd.x << " A" << std::endl;
   std::cout << "\ty = " << mat.ucd.y << " A" << std::endl;
   std::cout << "\tz = " << mat.ucd.z << " A" << std::endl;

   std::cout << std::endl;

   return EXIT_SUCCESS;

}

int read_coordinatefile() {

   std::string filename;

   if (mat.name == "ndfeti12")
      filename = "./coordinates/ndfeti12/" + std::to_string(sys.ti) + ".coords";
   else if (mat.name == "interface")
      filename = "./coordinates/interfaces/mirror_shifted.coords";
   else
      filename = "./coordinates/" + mat.name + ".coords";

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
      if (mat.name != "interface") tmp.mat = determine_element_id(tmp.element);
      else tmp.mat = determine_element_id_for_interface_system(tmp.element, tmp.pos.z);

      /* scale atom coordinates (except for interface system as they're already scaled) */
      if (mat.name != "interface") tmp.pos = tmp.pos * mat.ucd;

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

int modify_element_ids_for_interface_system() {



   for (int i=0; i<unitcell.size(); ++i) {

      /* this is the alpha fe */
      if (unitcell[i].pos.z > 74.0 &&
          unitcell[i].pos.z < 182.0 &&
          unitcell[i].mat == 1)

         unitcell[i].mat = 4;

      /* this is the second ndfeb phase */
      else if (unitcell[i].pos.z > 182.0) unitcell[i].mat += 5;

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

int neighbour_routines() {

   if (sys.rt_shell_cut) identify_re_neighbours();
   identify_tm_neighbours();

   return EXIT_SUCCESS;

}

/* this function identifies the atoms nearest to each rare-earth */
int identify_re_neighbours() {

   int start;
   int end;

   start = (supercell.size()-unitcell.size())/2;
   end   = (supercell.size()+unitcell.size())/2;

   double shell_1_distance = 1e3;
   double shell_2_distance = 1e3;
   double shell_3_distance = 1e3;
   double shell_4_distance = 1e3;
   double shell_5_distance = 1e3;

   std::ofstream re_neighbours ("re_neighbours.dat");

   std::cout << std::endl;

   /* loop through central cell of supercell */
   for (int i=start; i<end; ++i) {

      /* only continue with next loop if first is rare-earth */
      if (supercell[i].is_re())

         for (int j=0; j<supercell.size(); ++j) {

            /* define pair or atoms */
            pair_t pair;
            pair.i = supercell[i];
            pair.j = supercell[j];

            /* only interested in transition metal neighbours */
            if (pair.j.is_tm()) {

               double separation = pair.rij();

               re_neighbours << separation << std::endl;

               if (shell_1_distance-separation > 1e-10) {
                  shell_5_distance = shell_4_distance;
                  shell_4_distance = shell_3_distance;
                  shell_3_distance = shell_2_distance;
                  shell_2_distance = shell_1_distance;
                  shell_1_distance = separation;
               }

               else if (shell_2_distance-separation > 1e-10 && separation-shell_1_distance > 1e-10) {
                  shell_5_distance = shell_4_distance;
                  shell_4_distance = shell_3_distance;
                  shell_3_distance = shell_2_distance;
                  shell_2_distance = separation;
               }

               else if (shell_3_distance-separation > 1e-10 && separation-shell_2_distance > 1e-10) {
                  shell_5_distance = shell_4_distance;
                  shell_4_distance = shell_3_distance;
                  shell_3_distance = separation;
               }

               else if (shell_4_distance-separation > 1e-10 && separation-shell_3_distance > 1e-10) {
                  shell_5_distance = shell_4_distance;
                  shell_4_distance = separation;
               }

               else if (shell_5_distance-separation > 1e-10 && separation-shell_4_distance > 1e-10) {
                  shell_5_distance = separation;
               }

            }
         }
   }

   std::cout << "RE--TM neighbours:\n\n";
   std::cout << "\t1st neighbour distance: " << shell_1_distance << " A" << std::endl;
   std::cout << "\t2nd neighbour distance: " << shell_2_distance << " A" << std::endl;
   std::cout << "\t3rd neighbour distance: " << shell_3_distance << " A" << std::endl;
   std::cout << "\t4th neighbour distance: " << shell_4_distance << " A" << std::endl;
   std::cout << "\t5th neighbour distance: " << shell_5_distance << " A" << std::endl;
   std::cout << std::endl;

   if (sys.rt_shell == 1) sys.rcut_rt = (shell_2_distance + shell_1_distance)/2.0;
   else if (sys.rt_shell == 2) sys.rcut_rt = (shell_3_distance + shell_2_distance)/2.0;
   else if (sys.rt_shell == 3) sys.rcut_rt = (shell_4_distance + shell_3_distance)/2.0;
   else if (sys.rt_shell == 4) sys.rcut_rt = (shell_5_distance + shell_4_distance)/2.0;
   else if (sys.rt_shell == 5) sys.rcut_rt = shell_5_distance + 0.001;

   std::cout << "\tRE-TM cut-off set to " << sys.rcut_rt << " A\n";

   return EXIT_SUCCESS;

}

int identify_tm_neighbours() {

   int start;
   int end;

   start = (supercell.size()-unitcell.size())/2;
   end   = (supercell.size()+unitcell.size())/2;

   /* cut-off radii */
   double shell1cut = 1.0;
   double shell2cut = 2.0;
   double shell3cut = 3.0;
   double shell4cut = 4.0;
   double shell5cut = 5.0;

   /* atoms in each shell */
   double n_shell1_total = 0;
   double n_shell2_total = 0;
   double n_shell3_total = 0;
   double n_shell4_total = 0;
   double n_shell5_total = 0;

   double n_tm_atoms = 0;

   for (int i=start; i<end; ++i) {
      if (supercell[i].is_tm()) {

         n_tm_atoms ++;

         double n_shell1 = 0;
         double n_shell2 = 0;
         double n_shell3 = 0;
         double n_shell4 = 0;
         double n_shell5 = 0;

         for (int j=0; j<supercell.size(); ++j) {

            /* create a 'pair' */
            pair_t pair;
            pair.i = supercell[i];
            pair.j = supercell[j];

            if (pair.j.is_tm() && !pair.are_same_atom()) {

               if (pair.rij() <= shell1cut) {
                  n_shell1 ++;
                  n_shell2 ++;
                  n_shell3 ++;
                  n_shell4 ++;
                  n_shell5 ++;
               }

               else if (pair.rij() <= shell2cut) {
                  n_shell2 ++;
                  n_shell3 ++;
                  n_shell4 ++;
                  n_shell5 ++;
               }

               else if (pair.rij() <= shell3cut) {
                  n_shell3 ++;
                  n_shell4 ++;
                  n_shell5 ++;
               }

               else if (pair.rij() <= shell4cut) {
                  n_shell4 ++;
                  n_shell5 ++;
               }

               else if (pair.rij() <= shell5cut)
                  n_shell5 ++;
            }
         }

         n_shell1_total += n_shell1;
         n_shell2_total += n_shell2;
         n_shell3_total += n_shell3;
         n_shell4_total += n_shell4;
         n_shell5_total += n_shell5;

      }
   }

   std::cout << "\nTM--TM neighbours:\n\n";
   std::cout << "\t1.0 A shell: " << n_shell1_total / n_tm_atoms << std::endl;
   std::cout << "\t2.0 A shell: " << n_shell2_total / n_tm_atoms << std::endl;
   std::cout << "\t3.0 A shell: " << n_shell3_total / n_tm_atoms << std::endl;
   std::cout << "\t4.0 A shell: " << n_shell4_total / n_tm_atoms << std::endl;
   std::cout << "\t5.0 A shell: " << n_shell5_total / n_tm_atoms << std::endl;

   return EXIT_SUCCESS;

}
