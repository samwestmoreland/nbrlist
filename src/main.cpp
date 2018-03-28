#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <chrono>

/* local header files */
#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"
#include "./main.hpp"
#include "./initialise.hpp"

int main (int argc, char *argv[]) {

   auto now = std::chrono::system_clock::now();
   std::time_t current_time = std::chrono::system_clock::to_time_t(now);

   std::cout << std::endl;
   std::cout << "***************************************************\n\n";
   std::cout << "                   Nbrlist\n\n";
   std::cout << "       Compiled on " << __DATE__ << " at " << __TIME__ << "\n";
   std::cout << "       Run on " << std::ctime(&current_time) << "\n";
   std::cout << "***************************************************\n";

   /* global variables */
   int system_dimension = 0;
   int n_tracked_cells = 0;

   /* default parameters */
   sys.domainwall = false;
   sys.zrdoping = false;
   sys.zrconcentration = 0.0;

   /* read parameters from input file */
   parse_input("ucf_inputfile");

   std::cout << "\nparameters read from input file:\n";
   std::cout << "\tmaterial: " << sys.material << std::endl;
   std::cout << "\tcut-off radius: " << sys.rcut << std::endl;
   std::cout << "\tzr concentration: " << sys.zrconcentration << std::endl;

   // if (sys.tracking) {
   //    sys_dimension = ask_for_sys_dimensions(sys_dimension);
   //    n_tracked_cells = ask_for_n_tracked_cells(n_tracked_cells);
   // }

   /* exit program if material not recognised */
   if (sys.material_int > 7 || sys.material_int < 1) {
      std::cout << "invalid material. exiting.\n";
      exit(EXIT_SUCCESS);
   }

   else initialise_material(sys.material_int, sys.material);

   /* ************************
    * *** print parameters ***
    * ************************/

   // if (tracking)
   // {
   //     std::cout << "tracking system dimension: " << system_dimension << std::endl;
   //     std::cout << "no. of tracked cells: " << n_tracked_cells << std::endl;
   // }

   std::cout << std::endl;

   array_to_rasmol(unitcell, "unitcell");

//   /* if zr concentration != 0, the unitcell must be expanded to get even distribution */
//   if (system.zrconcentration != 0)
//      expand_unitcell_and_substitute_zr_atoms(system.zrconcentration);

   /* this function generates a supercell, and using that populates an interactions array */
   calculate_interactions();

   /**************************************/
   /*** calculate species interactions ***/
   /**************************************/

   /* vector to hold interaction energy of each species with every other species */
   std::vector<std::vector<double> > species_interactions;
   species_interactions.resize(materials.size());
   for (int i=0; i<materials.size(); ++i) species_interactions[i].resize(materials.size());

   /* vector with same dimensions to hold number of interactions */
   std::vector<std::vector<int> > n_interactions;
   n_interactions.resize(materials.size());
   for (int i=0; i<materials.size(); ++i) n_interactions[i].resize(materials.size());

   /* initialise both arrays to zero */
   for (int i=0; i<materials.size(); ++i)
       for (int j=0; j<materials.size(); ++j) {

           species_interactions[i][j] = 0;
           n_interactions[i][j] = 0;

       }

   for (int i=0; i<uc_interactions.size(); ++i) {

       int j = uc_interactions[i].j.mat;
       int k = uc_interactions[i].i.mat;

       species_interactions[j][k] += uc_interactions[i].exchange;
       n_interactions[j][k] ++;
   }

   /**************************************/
   /*********** output arrays ************/
   /**************************************/

   std::cout << "\nn_interactions matrix\n";

   for (int i=0; i<materials.size(); ++i) {
       for (int j=0; j<materials.size(); ++j) {

           std::cout << n_interactions[i][j] << "\t";

       }

       std::cout << std::endl;
   }

   std::cout << "\nexchange matrix (mean)\n";

   for (int i=0; i<materials.size(); ++i) {
       for (int j=0; j<materials.size(); ++j) {
          if (n_interactions[i][j] != 0) {

             species_interactions[i][j] /= n_interactions[i][j];
             std::cout << species_interactions[i][j] << "\t";

          }

          else std::cout << "0.000000000" << "\t";
       }

       std::cout << std::endl;
   }

   std::cout << std::endl;

   /* generate a domain wall system */
   if (sys.domainwall == true) generate_domain_wall_system(sys.dw_dim);

   /***********************************************/
   /*** generate large system for cell tracking ***/
   /***********************************************/

   // if (tracking) {

   //    /* array to hold atoms in large system */
   //    // std::vector < std::vector < std::vector < std::vector <atom_t> > > > system;

   //    // generate_large_system(uc_interactions, unitcell, system, ucd, n_tracked_cells, materials.size(), system_dimension);

   // }

   return 0;

}

int convert_material_string_to_integer(std::string const& material) {

   int material_int = 0;

   if (material == "bccfe") material_int = 1;
   else if (material == "ndfeb") material_int = 2;
   else if (material == "ndfe12") material_int = 3;
   else if (material == "smfe12") material_int = 4;
   else if (material == "smzrfe12") material_int = 5;
   else if (material == "interface") material_int = 6;
   else if (material == "interface_mirror") material_int = 7;

   return material_int;
}

// for now this function will include tracked cell calculation */
int generate_large_system(std::vector<int_t>& uc_interactions,
                          std::vector<atom_t>& unitcell,
                          std::vector < std::vector < std::vector < std::vector <atom_t> > > >& system,
                          vec_t ucd,
                          int n_tracked_cells,
                          int n_materials,
                          int system_dimension) {

   /* system dimensions in unit cells */
   vec_t sd;
   sd.x = system_dimension;
   sd.y = system_dimension;
   sd.z = system_dimension;

   std::cout << "\ndimensions of large system [A]: ("
   << sd.x * ucd.x << ", "
   << sd.y * ucd.y << ", "
   << sd.z * ucd.z << ")\n";

   int ns = sd.x * sd.y * sd.z * unitcell.size(); // number of atoms in system
   std::cout << "number of atoms in system: " << ns << "\n";

   std::ofstream largeucf ("large.ucf");

   /* output coordinates to unit cell file */
   largeucf << "# Unit cell size:\n"
            <<  sd.x*ucd.x << "\t"
            <<  sd.y*ucd.y << "\t"
            <<  sd.z*ucd.z << "\n"
            << "# Unit cell vectors:\n"
            << "1.0  0.0  0.0\n"
            << "0.0  1.0  0.0\n"
            << "0.0  0.0  1.0\n"
            << "# Atoms num, id cx cy cz mat lc hc\n"
            << ns << "\n";

   /* global id counter */
   int gid_counter = 0;
   int tracked = 0;

   /* open file for rasmol output of system */
   std::ofstream sysmol ("system.xyz");
   sysmol << ns << "\n\n";

   /* resize vectors */
   system.resize(sd.x);
   for (int i=0; i<sd.x; i++) {
      system[i].resize(sd.y);
      for (int j=0; j<sd.y; j++) {
         system[i][j].resize(sd.z);
         for (int k=0; k<sd.z; k++) {

            /* loop through unitcell atoms */
            for (int atom=0; atom<unitcell.size(); ++atom) {

               atom_t temp;
               vec_t uc;
               uc.x = i;
               uc.y = j;
               uc.z = k;

               temp.aid = unitcell[atom].aid;
               temp.gid = gid_counter;
               gid_counter ++;

               temp.element = unitcell[atom].element;

               temp.mat = unitcell[atom].mat;

               /*** cell tracking ***/

               if (  (i <  (sd.x+n_tracked_cells)*0.5)
                  && (i >= (sd.x-n_tracked_cells)*0.5)

                  && (j <  (sd.y+n_tracked_cells)*0.5)
                  && (j >= (sd.y-n_tracked_cells)*0.5)

                  && (k <  (sd.z+n_tracked_cells)*0.5)
                  && (k >= (sd.z-n_tracked_cells)*0.5) ) {

                  temp.mat += n_materials;
                  temp.element = 'h';
                  tracked ++;
               }

               temp.pos = unitcell[atom].pos + uc*ucd;

               /* output coordinates to file */
               sysmol << temp.element << "\t"
                      << temp.pos.x << "\t"
                      << temp.pos.y << "\t"
                      << temp.pos.z << "\n";

               /* calculate atom coordinates within large system */
               vec_t sys_co;
               sys_co.x = temp.pos.x / double(ucd.x) / double(sd.x);
               sys_co.y = temp.pos.y / double(ucd.y) / double(sd.y);
               sys_co.z = temp.pos.z / double(ucd.z) / double(sd.z);

               /* output to unit cell file */
               largeucf
                  << temp.gid << "\t"
                  << sys_co.x << "\t"
                  << sys_co.y << "\t"
                  << sys_co.z << "\t"
                  << temp.mat << "\t"
                  << 0 << "\t"
                  << 0 << "\n";

               system[i][j][k].push_back(temp);

            }
         }
      }
   }

   std::cout << "total tracked atoms: " << tracked << "\n";
   std::cout << "total tracked cells (3D): " << tracked/unitcell.size() << "\n";

   largeucf << "# interactions n exctype, id i j dx dy dz Jij\n";

   /*******************************************************/
   /**  determine interactions for every atom in system  **/
   /**         using pre-calculated neighbour list       **/
   /*******************************************************/

   /* initialise array to hold interactions for whole system */
   std::vector<int_t> interactions;

   /* interaction counter */
   int int_counter = 0;

   /* loop through every atom in system */
   for (int i=0; i<sd.x; i++)
   for (int j=0; j<sd.y; j++)
   for (int k=0; k<sd.z; k++)

   /* first atom */
   for (int atom=0; atom<unitcell.size(); atom++)

   /* loop through interaction information */
   for (int p=0; p<uc_interactions.size(); p++) {

      /* if interaction info refers to correct atom */
      if (system[i][j][k][atom].aid == uc_interactions[p].i.aid) {

         int_t temp;

         temp.iid = int_counter;
         temp.i.gid = system[i][j][k][atom].gid;

         temp.i.mat = system[i][j][k][atom].mat;
         temp.i.element = system[i][j][k][atom].element;

         temp.j.element = uc_interactions[p].j.element;

         temp.exchange = uc_interactions[p].exchange;

         /* assume atom j is within system to begin with */
         temp.disp.x = 0;
         temp.disp.y = 0;
         temp.disp.z = 0;

         /* check if atom j is within system boundaries */
         int ucx = i + uc_interactions[p].disp.x;
         int ucy = j + uc_interactions[p].disp.y;
         int ucz = k + uc_interactions[p].disp.z;

         /* if any of these conditions satisfied
          * then atom is out of bounds
          */

         if (ucx < 0 || ucx >= sd.x ||
             ucy < 0 || ucy >= sd.y ||
             ucz < 0 || ucz >= sd.z  ) {

            /* periodic boundaries conditions */
            if ( ucx < 0 ) {
               ucx += sd.x;
               temp.disp.x = -1;

                  //    // determine mat of atom j
                  //    if (temp.i.mat==5)
                  //    {
                  //        if (temp.j.element=="Nd") temp.j.mat = 5;
                  //        else temp.j.mat = 6;
                  //    }
                  //
                  //    else if (temp.i.mat==6)
                  //    {
                  //        if (temp.j.element=="Fe") temp.j.mat = 6;
                  //        else temp.j.mat = 5;
                  //    }

            }

            if ( ucy < 0 ) {
                 ucy += sd.y;
                 temp.disp.y = -1;
            }

            if ( ucz < 0 ) {
                 ucz += sd.z;
                 temp.disp.z = -1;
            }

            if ( ucx >= sd.x ) {
                 ucx -= sd.x;
                 temp.disp.x = 1;
            }

            if ( ucy >= sd.y ) {
                 ucy -= sd.y;
                 temp.disp.y = 1;
            }

            if ( ucz >= sd.z ) {
                 ucz -= sd.z;
                 temp.disp.z = 1;
            }

            /* having changed uc coordinates obtain j.gid */
            temp.j.gid = system[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

         }

         /* else it is within bounds so simply extract j.gid */
         else temp.j.gid = system[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

         interactions.push_back(temp);

         /* increment interaction id */
         int_counter ++;

      }

      std::cout << "number of interactions in system: " << interactions.size() << "\n\n";

      largeucf << interactions.size() << "\tisotropic\n";

      // output interaction info to file
      for (int i=0; i<interactions.size(); i++)

      largeucf << interactions[i].iid << "\t"
          << interactions[i].i.gid << "\t"
          << interactions[i].j.gid << "\t"
          << interactions[i].disp.x << "\t"
          << interactions[i].disp.y << "\t"
          << interactions[i].disp.z << "\t"
          << interactions[i].exchange << "\n";
   }

   largeucf.close();

   return EXIT_SUCCESS;
}

vec_t get_uc_dimensions_from_zr_content(double zrconcentration) {

    double atom_rad = (zrconcentration - 8.22624)/(-4.51952);
    double a = 0.0929088 * atom_rad + 0.830892;
    double c = -0.00593675 * zrconcentration + 1;

    a *= 8.497;
    c *= 4.687;

    vec_t ucd;
    ucd.x = a;
    ucd.y = a;
    ucd.z = c;
    return ucd;
}

/* returns a material */
int determine_material_id(std::string const& in_material, std::vector<material_t>& materials) {

   /* define out_material which is a material_t to pushback into material array */
   material_t out_material;
   out_material.name = in_material;
   out_material.id = 0;

   /* if this is the first time we've encountered this material */
   if (materials.size() == 0) {

      materials.push_back(out_material);
      return out_material.id;

   }

   for (int i=0; i<materials.size(); ++i) {

      if (in_material == materials[i].name) {

         out_material.id = materials[i].id;

         return out_material.id;

      }
   }

   /* if we get to here, material has not been found in material array */
   out_material.id = materials.size();
   materials.push_back(out_material);

   return out_material.id;
}

void populate_supercell() {

   /* loop through dimensions */
   for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
         for (int k=0; k<3; k++)

            /* loop through atoms in unitcell */
            for (int atom=0; atom<unitcell.size(); atom++) {

               atom_t temp;
               vec_t uc;
               uc.x = i;
               uc.y = j;
               uc.z = k;

               temp.aid = unitcell[atom].aid;
               temp.element = unitcell[atom].element;
               temp.mat = unitcell[atom].mat;

               /* replicate unitcell atoms */
               temp.pos = unitcell[atom].pos + (uc * ucd);

               /* label unitcell coordinates */
               temp.uc = uc;

               /* dummy values for unneeded struct elements */
               temp.hcat = 0;
               temp.gid = 0;   // this isn't needed yet as the calculations are generated for a large system later

               /* place atom in array */
               supercell.push_back(temp);
            }

   std::cout << "atoms in super cell: "
             << supercell.size() << std::endl;

   array_to_rasmol(supercell, "supercell");

}

int calculate_interactions() {

   std::cout << "\ncreating and populating super cell...\n";

   populate_supercell();

   std::cout << "calculating interactions...\n\n";

   if (sys.material_int == 2)
      std::cout <<
         "TM-TM exchange factor = " << sys.tt_factor << "\n" <<
         "RE-TM exchange factor = " << sys.rt_factor << "\n";

   int start;
   int end;

   int interaction_count = 0;
   int fefe_interaction_count = 0;
   int refe_interaction_count = 0;

   double total_neighbour_distance = 0;
   double total_fefe_neighbour_distance = 0;
   double total_refe_neighbour_distance = 0;


      start = (supercell.size()-unitcell.size())/2;
      end   = (supercell.size()+unitcell.size())/2;

      /* loop through central cell */
      for (int i=start; i<end; ++i) {

         /* loop through super cell (looking for atom j) */
         for (int j=0; j<supercell.size(); ++j) {

            /* calculate interatomic distance */
            double rij = calculate_rij(supercell[i].pos, supercell[j].pos);

            /* if distance less than rcut and not same atom */
            if (rij < sys.rcut && rij > 1e-30) {

               /* add neighbour distance to total */
               total_neighbour_distance += rij;

               /* create an interaction */
               int_t temp;
               temp.i = supercell[i];
               temp.j = supercell[j];

               temp.iid = interaction_count;

               /* unitcell displacement */
               temp.disp = temp.j.uc - temp.i.uc;

               /* calculate exchange energy */
               temp.exchange = calculate_jij(temp.i.element, temp.j.element, rij, sys.tt_factor, sys.rt_factor, sys.material_int);

               /* put interaction into array */
               if (temp.exchange!=0) {

                  uc_interactions.push_back(temp);
                  interaction_count ++;

                  if (
                        (temp.i.element == "Fe8i" || temp.i.element == "Fe8j" || temp.i.element == "Fe8f" || temp.i.element == "Fe") &&
                        (temp.j.element == "Fe8i" || temp.j.element == "Fe8j" || temp.j.element == "Fe8f" || temp.j.element == "Fe")
                     ) {

                     fefe_interaction_count ++;
                     total_fefe_neighbour_distance += rij;
                  }

                  else if (
                        ((temp.i.element == "Fe8i" || temp.i.element == "Fe8j" || temp.i.element == "Fe8f" || temp.i.element == "Fe") &&
                         (temp.j.element == "Sm" || temp.j.element == "Nd"))
                        ||
                        ((temp.j.element == "Fe8i" || temp.j.element == "Fe8j" || temp.j.element == "Fe8f" || temp.j.element == "Fe") &&
                         (temp.i.element == "Sm" || temp.i.element == "Nd"))
                        ) {
                     refe_interaction_count ++;
                     total_refe_neighbour_distance += rij;
                  }
               }

            }
         }
      }

      std::cout
         << "total interactions: "
         << uc_interactions.size() << std::endl;

      std::cout
         << "fe-fe interactions: "
         << fefe_interaction_count << std::endl;

      std::cout
         << "re-fe interactions: "
         << refe_interaction_count << std::endl;

      std::cout << std::endl;

      std::cout
         << "mean neighbour distance: "
         << total_neighbour_distance/interaction_count
         << " A" << std::endl;

      std::cout
         << "mean fe-fe neighbour distance: "
         << total_fefe_neighbour_distance/fefe_interaction_count
         << " A" << std::endl;

      std::cout
         << "mean re-fe neighbour distance: "
         << total_refe_neighbour_distance/refe_interaction_count
         << " A" << std::endl;

      std::cout << std::endl;

      outfile
         << "# Interactions n exctype, id i j dx dy dz Jij\n"
         << uc_interactions.size() << "\tisotropic\n";

      for (int i=0; i<uc_interactions.size(); ++i)
         outfile
            << uc_interactions[i].iid       << "\t"
            << uc_interactions[i].i.aid     << "\t"
            << uc_interactions[i].j.aid     << "\t"
            << uc_interactions[i].disp.x    << "\t"
            << uc_interactions[i].disp.y    << "\t"
            << uc_interactions[i].disp.z    << "\t"
            << uc_interactions[i].exchange  << "\n";

      std::cout << "interaction data output to \'output.ucf\'\n";

      return EXIT_SUCCESS;
}

/* calculate exchange energy */
double calculate_jij(std::string const& i_type,
                     std::string const& j_type,
                     double rij,
                     double tt_factor,
                     double rt_factor,
                     int exchange_fn) {

   switch (exchange_fn) {

      case 1 : { // bccfe
         /* parameters from Pajda (2001) */
         /* factor 2 to correct for Hamiltonian */
         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;
         return c*(a/(rij*rij*rij)-b);
      }
         break;

      case 2 : // ndfeb

         return jij_ndfeb(i_type, j_type, rij);
         break;

      case 3 : { // ndfe12

         const double Fe_ratio = 1.2;      // loop from 1-2
         const double J0Nd = Fe_ratio * 4.06835e-20/16.0;

         /* R-Fe exchange */
         double RFeFraction = 0.2;         // loop from 0.0-1.0

         /* Fe-Fe (cutoff at r = 5.74A) */
         if (i_type == "Fe" && j_type == "Fe") {

            if ( rij <= 5.0 ) {

               double a = 36.9434;
               double b = 1.25094;
               double c = -0.229572;

               return 2.0*2.179872e-21 * (a*exp(-b*rij)+c) * Fe_ratio;
            }

            else return 0.0;
         }

         /* Nd-Fe (step function at r = 4A) */
         else if ((i_type=="Fe" && j_type=="Nd") || (i_type=="Nd" && j_type=="Fe")) {
            if (rij <= 4.0) return RFeFraction * J0Nd;       // ndfeb

            else return 0.0;
      	}

    		/* Nd-Nd */
      	else if ( i_type == "Nd" && j_type == "Nd" ) return 0.0;

      	else return 0.0;
    	}
      	break;

      case 4 : { /* smfe12 */

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         double smfe12_tt_factor = 0.65;
         double smfe12_rt_factor = 0.20;

         /* Sm-Sm */
         if(i_type=="Sm" && j_type=="Sm") return 0.0;

         /* Sm-Fe */
         if ((i_type=="Fe8i" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8i") ||
             (i_type=="Fe8j" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8j") ||
             (i_type=="Fe8f" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8f") ) {

            if (rij<=4.0) return smfe12_rt_factor * c*(a/(rij*rij*rij)-b);
            else return 0.0;
         }

         // Fe-Fe (cutoff at r = 5.74A)
         else if ((i_type=="Fe8i" && j_type=="Fe8i") ||
                  (i_type=="Fe8i" && j_type=="Fe8j") ||
                  (i_type=="Fe8i" && j_type=="Fe8f") ||
                  (i_type=="Fe8j" && j_type=="Fe8i") ||
                  (i_type=="Fe8j" && j_type=="Fe8j") ||
                  (i_type=="Fe8j" && j_type=="Fe8f") ||
                  (i_type=="Fe8f" && j_type=="Fe8i") ||
                  (i_type=="Fe8f" && j_type=="Fe8j") ||
                  (i_type=="Fe8f" && j_type=="Fe8f"))

            return smfe12_tt_factor * c*(a/(rij*rij*rij)-b);

         else return 0.0;
      }
         break;

      case 5 : { /* smzrfe12 */

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         /* Sm-Sm */
         if(i_type=="Sm" && j_type=="Sm") return 0.0;

         /* Sm-Fe */
         if ((i_type=="Fe8i" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8i") ||
             (i_type=="Fe8j" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8j") ||
             (i_type=="Fe8f" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8f") ) {

            if (rij<=4.0) return rt_factor * c*(a/(rij*rij*rij)-b);
            else return 0.0;
         }

         // Fe-Fe (cutoff at r = 5.74A)
         else if ((i_type=="Fe8i" && j_type=="Fe8i") ||
                  (i_type=="Fe8i" && j_type=="Fe8j") ||
                  (i_type=="Fe8i" && j_type=="Fe8f") ||
                  (i_type=="Fe8j" && j_type=="Fe8i") ||
                  (i_type=="Fe8j" && j_type=="Fe8j") ||
                  (i_type=="Fe8j" && j_type=="Fe8f") ||
                  (i_type=="Fe8f" && j_type=="Fe8i") ||
                  (i_type=="Fe8f" && j_type=="Fe8j") ||
                  (i_type=="Fe8f" && j_type=="Fe8f"))

            return tt_factor * c*(a/(rij*rij*rij)-b);

         else return 0.0;
      }
         break;

      case 6 : /* interface */

         return jij_ndfeb(i_type, j_type, rij);
         break;

      case 7 : /* interface mirror */

         return jij_ndfeb(i_type, j_type, rij);
         break;

      default :
         std::cout << "invalid option. exiting.\n";
         exit(EXIT_FAILURE);

   } /* end of switch */

} /* end of calculate_jij function */

double jij_ndfeb(std::string const& i_type,
                 std::string const& j_type,
                 double rij) {

   /* richard's code */
   const double a = 36.9434;
   const double b = 1.25094;
   const double c = -0.229572;
   const double Fe_ratio_ndfeb = 0.69618016759 * 1.07692307692; // 560/520 = 1.07692307692
   const double J0Nd_ndfeb = Fe_ratio_ndfeb * 4.06835e-20 / 16.0;

   /* nd-nd */
   if (i_type == "Nd" && j_type == "Nd") return 0.0;

   /* nd-fe (step function at r = 4A) */
   else if ((i_type == "Fe" && j_type == "Nd") || (i_type == "Nd" && j_type == "Fe")) {
      if (rij <= 4.0) return 0.33 * J0Nd_ndfeb;       // ndfeb
      else return 0.0;
   }

   /* fe-fe (cutoff at r = 5.74A) */
   else if (i_type == "Fe" && j_type == "Fe") {
      // if(rij<=5.0) return -2.0*2.179872e-21*(A*exp(-B*rij)+C);
      // correct for Tc = 600
      if (rij <= 5.0) return 2.0 * 2.179872e-21 * (a*exp(-b*rij)+c) * Fe_ratio_ndfeb;
      else return 0.0;
   }

   /* boron */
   else return 0.0;
}

/* function to calculate distance */
double calculate_rij (vec_t& i, vec_t& j) {

    vec_t d = j-i;
    double distance = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);

    return distance;
}

int generate_domain_wall_system(vec_t dw_dim) {

   std::cout << "initialising domain wall system\n\n";

   /* determine dimensions of system in unitcells */
   /* for long system in z
      x = 2 nm
      y = 2 nm
      z = 12 nm */

   std::ofstream dwucf ("domainwall.ucf");
   vec_t sd; /* system dimensions in unitcells (dw_dim comes from input file and is in nm) */


   sd.x = floor(dw_dim.x*10.0/ucd.x+0.5);
   sd.y = floor(dw_dim.y*10.0/ucd.y+0.5);
   sd.z = floor(dw_dim.z*10.0/ucd.z+0.5);

   std::cout << "dimensions of domain wall system: ";
   std::cout << dw_dim.x << " x " << dw_dim.y << " x " << dw_dim.z << " nm\n";

   int natoms = sd.x * sd.y * sd.z * unitcell.size(); // number of atoms in system
   std::cout << "number of atoms in system: " << natoms << "\n";

   /* output coordinates to unit cell file */
   dwucf << "# Unit cell size:\n"
      <<  sd.x*ucd.x << "\t"
      <<  sd.y*ucd.y << "\t"
      <<  sd.z*ucd.z << "\n"
      << "# Unit cell vectors:\n"
      << "1.0  0.0  0.0\n"
      << "0.0  1.0  0.0\n"
      << "0.0  0.0  1.0\n"
      << "# Atoms num, id cx cy cz mat lc hc\n"
      << natoms << "\n";

   /* global id counter */
   int gid_counter = 0;

   /* open file for rasmol output of system */
   std::ofstream sysmol ("domainwall.xyz");
   sysmol << natoms << "\n\n";

   /* resize vectors */
   domainwallsystem.resize(sd.x);
   for (int i=0; i<sd.x; i++) {
      domainwallsystem[i].resize(sd.y);
      for (int j=0; j<sd.y; j++) {
         domainwallsystem[i][j].resize(sd.z);
         for (int k=0; k<sd.z; k++) {

            /* loop through unitcell atoms */
            for (int atom=0; atom<unitcell.size(); ++atom) {

               atom_t tmp;
               vec_t uc;
               uc.x = i;
               uc.y = j;
               uc.z = k;

               tmp.aid = unitcell[atom].aid;
               tmp.gid = gid_counter;

               tmp.element = unitcell[atom].element;

               tmp.mat = unitcell[atom].mat;

               tmp.pos = unitcell[atom].pos + uc*ucd;
               gid_counter ++;

               tmp.hcat = i;

               /* calculate atom coordinates within large system */
               vec_t sys_coord;
               sys_coord.x = tmp.pos.x / double(ucd.x) / double(sd.x);
               sys_coord.y = tmp.pos.y / double(ucd.y) / double(sd.y);
               sys_coord.z = tmp.pos.z / double(ucd.z) / double(sd.z);

               /* we want to split the system into two halves */
               if (uc.x > (sd.x-1)/2) {
                  tmp.mat += materials.size();
               }

               /* output to unit cell file */
               dwucf << tmp.gid << "\t"
                  << sys_coord.x << "\t"
                  << sys_coord.y << "\t"
                  << sys_coord.z << "\t"
                  << tmp.mat << "\t"
                  << 0 << "\t"
                  << tmp.hcat << "\n";

               domainwallsystem[i][j][k].push_back(tmp);
            }
         }
      }
   }

   dwucf << "# interactions n exctype, id i j dx dy dz Jij\n";

   /*******************************************************/
   /**  determine interactions for every atom in system  **/
   /**         using pre-calculated neighbour list       **/
   /*******************************************************/

   /* initialise array to hold interactions for whole system */
   std::vector<int_t> interactions;

   /* interaction counter */
   int int_counter = 0;

   int unitcell_size;
   unitcell_size = unitcell.size();

   /* loop through every atom in system */
   for (int i=0; i<sd.x; i++) {
      for (int j=0; j<sd.y; j++) {
         for (int k=0; k<sd.z; k++) {

            /* first atom */
            for (int atom=0; atom<unitcell_size; atom++) {

               /* loop through interaction information */
               for (int p=0; p<uc_interactions.size(); p++) {

                  /* if interaction info refers to correct atom */
                  if (domainwallsystem[i][j][k][atom].aid == uc_interactions[p].i.aid) {

                     int_t tmp;

                     tmp.iid = int_counter;
                     tmp.i.gid = domainwallsystem[i][j][k][atom].gid;

                     tmp.i.mat = domainwallsystem[i][j][k][atom].mat;
                     tmp.i.element = domainwallsystem[i][j][k][atom].element;

                     tmp.j.element = uc_interactions[p].j.element;

                     tmp.exchange = uc_interactions[p].exchange;

                     /* assume atom j is within system to begin with */
                     tmp.disp.x = 0;
                     tmp.disp.y = 0;
                     tmp.disp.z = 0;

                     /* check if atom j is within system boundaries */
                     int ucx = i + uc_interactions[p].disp.x;
                     int ucy = j + uc_interactions[p].disp.y;
                     int ucz = k + uc_interactions[p].disp.z;

                     /* if any of these conditions satisfied
                      * then atom is out of bounds */

                     if (ucx < 0 || ucx >= sd.x ||
                         ucy < 0 || ucy >= sd.y ||
                         ucz < 0 || ucz >= sd.z ) {

                        /* periodic boundaries conditions */
                        if ( ucx < 0 ) {
                           tmp.exchange *= -1; /* antiferro periodic boundaries */
                           ucx += sd.x;
                           tmp.disp.x = -1;

                           //    // determine mat of atom j
                           //    if (tmp.i.mat==5)
                           //    {
                           //        if (tmp.j.element=="Nd") tmp.j.mat = 5;
                           //        else tmp.j.mat = 6;
                           //    }
                           //
                           //    else if (tmp.i.mat==6)
                           //    {
                           //        if (tmp.j.element=="Fe") tmp.j.mat = 6;
                           //        else tmp.j.mat = 5;
                           //    }

                        }

                        if ( ucy < 0 ) {
                           ucy += sd.y;
                           tmp.disp.y = -1;
                        }

                        if ( ucz < 0 ) {
                           ucz += sd.z;
                           tmp.disp.z = -1;
                        }

                        if ( ucx >= sd.x ) {
                           tmp.exchange *= -1; /* antiferro periodic boundaries */
                           ucx -= sd.x;
                           tmp.disp.x = 1;
                        }

                        if ( ucy >= sd.y ) {
                           ucy -= sd.y;
                           tmp.disp.y = 1;
                        }

                        if ( ucz >= sd.z ) {
                           ucz -= sd.z;
                           tmp.disp.z = 1;
                        }

                        /* having changed uc coordinates obtain j.gid */
                        tmp.j.gid = domainwallsystem[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

                     }

                     /* else it is within bounds so simply extract j.gid */
                     else tmp.j.gid = domainwallsystem[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

                     interactions.push_back(tmp);

                     /* increment interaction id */
                     int_counter ++;

                  }
               }
            }
         }
      }
   }

   // std::cout << "number of interactions in system: " << interactions.size() << "\n\n";

   dwucf << interactions.size() << "\tisotropic\n";

   // output interaction info to file
   for (int i=0; i<interactions.size(); i++)

      dwucf << interactions[i].iid << "\t"
            << interactions[i].i.gid << "\t"
            << interactions[i].j.gid << "\t"
            << interactions[i].disp.x << "\t"
            << interactions[i].disp.y << "\t"
            << interactions[i].disp.z << "\t"
            << interactions[i].exchange << "\n";

   std::cout << "number of interactions found: " << interactions.size() << std::endl;
   std::cout << "interaction data output to file 'domainwall.ucf'\n";

   dwucf.close();

   return EXIT_SUCCESS;
}
