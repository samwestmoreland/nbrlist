#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <chrono>

/* project header files */
#include "./classes.hpp"
#include "./data.hpp"
#include "./io.hpp"
#include "./main.hpp"
#include "./initialise.hpp"
#include "./interactions.hpp"

int main (int argc, char *argv[]) {

   output_header();

   set_default_parameters();

   read_inputfile();

   get_unitcell_dimensions();

   read_coordinatefile();

   if (mat.name == "interface") populate_supercell_2D();
   else populate_supercell();

   if (mat.name != "interface") neighbour_routines();

   /* this function generates a supercell, and using that populates an interactions array */
   calculate_interactions();

   /* generate a domain wall system */
   if (sys.domainwall == true) generate_domain_wall_system();

   return 0;

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
int determine_element_id(std::string const& in_element) {

   /* element id is object to be returned */
   int element_id = 0;

   /* if no element has yet been found */
   if (elements.size() == 0) {

      elements.push_back(in_element);
      element_specific_atom_count.push_back(1);
      return element_id;

   }

   /* if we get to here, elements have been found previously */
   for (int i=0; i<elements.size(); ++i) {

      if (in_element == elements[i]) {
         element_id = i;
         element_specific_atom_count[i] ++;

         return element_id;
      }
   }

   /* if we get to here, element has not been found so add to array */
   elements.push_back(in_element);
   element_specific_atom_count.push_back(1);

   return elements.size()-1;
}

int determine_element_id_for_interface_system(std::string const& in_element, double zpos) {

   elements.resize(8);
   element_specific_atom_count.resize(8);

   elements[0] = "B";
   elements[1] = "Fe";
   elements[2] = "Nd";
   elements[3] = "Fe";
   elements[4] = "Co";
   elements[5] = "B";
   elements[6] = "Fe";
   elements[7] = "Nd";

   if (zpos < 74.0) {

      if (in_element == "B") {

         element_specific_atom_count[0] ++;
         return 0;

      }

      else if (in_element == "Fe") {

         element_specific_atom_count[1] ++;
         return 1;

      }

      else if (in_element == "Nd") {

         element_specific_atom_count[2] ++;
         return 2;

      }

      else {

         std::cout
            << "\nfatal error: foreign element: \""
            << in_element
            << "\" with z coordinate "
            << zpos
            << "\" found in first ndfeb phase. exiting.\n";

         exit(EXIT_FAILURE);

      }

   }

   else if (zpos >= 74.0 && zpos <= 182.0) {

      if (in_element == "Fe") {

         element_specific_atom_count[3] ++;
         return 3;

      }

      else if (in_element == "Co") {

         element_specific_atom_count[4] ++;
         return 4;

      }

      else {

         std::cout
            << "\nfatal error: foreign element: \""
            << in_element
            << "\" with z coordinate "
            << zpos
            << " found in alpha fe phase. exiting.\n";

         exit(EXIT_FAILURE);

      }

   }

   else if (zpos > 182.0) {

      if (in_element == "B") {

         element_specific_atom_count[5] ++;
         return 5;

      }

      else if (in_element == "Fe") {

         element_specific_atom_count[6] ++;
         return 6;

      }

      else if (in_element == "Nd") {

         element_specific_atom_count[7] ++;
         return 7;

      }

      else {

         std::cout
            << "\nfatal error: foreign element: \""
            << in_element
            << "\" with z coordinate "
            << zpos
            << "\" found in second ndfeb phase. exiting.\n";

         exit(EXIT_FAILURE);

      }

   }

   else {

      std::cout
         << "fatal error: found element with z coordinate "
         << zpos
         << " which is outside the correct range. exiting.\n";

      exit(EXIT_FAILURE);

   }

   return EXIT_SUCCESS;

}

/* function to calculate distance */
double calculate_rij (vec_t& i, vec_t& j) {

    vec_t d = j-i;
    double distance = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);

    return distance;
}

int generate_domain_wall_system() {

   std::cout << "initialising domain wall system\n\n";

   /* determine dimensions of system in unitcells */
   /* for long system in z
      x = 2 nm
      y = 2 nm
      z = 12 nm */

   std::ofstream dwucf ("domainwall.ucf");

   std::cout << "dimensions of domain wall system: ";
   std::cout
      << sys.dw_dim.x << " x "
      << sys.dw_dim.y << " x "
      << sys.dw_dim.z << " nm\n";

   /* system dimensions in unitcells (dw_dim comes from input file and is in nm) */
   vec_t sd;
   sd.x = floor(sys.dw_dim.x*10.0/mat.ucd.x+0.5);
   sd.y = floor(sys.dw_dim.y*10.0/mat.ucd.y+0.5);
   sd.z = floor(sys.dw_dim.z*10.0/mat.ucd.z+0.5);

   int natoms = sd.x * sd.y * sd.z * unitcell.size(); // number of atoms in system
   std::cout << "number of atoms in system: " << natoms << "\n";

   /* output coordinates to unit cell file */
   dwucf
      << "# Unit cell size:\n"
      <<  sd.x*mat.ucd.x << "\t"
      <<  sd.y*mat.ucd.y << "\t"
      <<  sd.z*mat.ucd.z << "\n"
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

   if (sys.centrepin == false) {

   /* resize vectors */
   dwsystem.resize(sd.x);
   for (int i=0; i<sd.x; i++) {
      dwsystem[i].resize(sd.y);
      for (int j=0; j<sd.y; j++) {
         dwsystem[i][j].resize(sd.z);
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

               tmp.pos = unitcell[atom].pos + uc * mat.ucd;
               gid_counter ++;

               tmp.hcat = i;

               /* calculate atom coordinates within large system */
               vec_t sys_coord;
               sys_coord.x = tmp.pos.x / double(mat.ucd.x) / double(sd.x);
               sys_coord.y = tmp.pos.y / double(mat.ucd.y) / double(sd.y);
               sys_coord.z = tmp.pos.z / double(mat.ucd.z) / double(sd.z);

               /* we want to split the system into two halves */
               if (uc.z > (sd.z-1)/2) {
                  tmp.mat += elements.size();
               }

               /* output to unit cell file */
               dwucf << tmp.gid << "\t"
                  << sys_coord.x << "\t"
                  << sys_coord.y << "\t"
                  << sys_coord.z << "\t"
                  << tmp.mat << "\t"
                  << 0 << "\t"
                  << tmp.hcat << "\n";

               if (tmp.mat == 0) sysmol << "H";
               else sysmol << "Ag";

               sysmol << "\t" << tmp.pos.x
                      << "\t" << tmp.pos.y
                      << "\t" << tmp.pos.z << "\n";

               dwsystem[i][j][k].push_back(tmp);
            }
         }
      }
   }
   }

   else if (sys.centrepin == true) {

   std::ofstream pinxyz ("pin.xyz");

   int n_materials = elements.size();
   n_materials *= 2;
   int n_pin = 0;

   /* resize vectors */
   dwsystem.resize(sd.x);
   for (int i=0; i<sd.x; i++) {
      dwsystem[i].resize(sd.y);
      for (int j=0; j<sd.y; j++) {
         dwsystem[i][j].resize(sd.z);
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

               tmp.pos = unitcell[atom].pos + uc * mat.ucd;
               gid_counter ++;

               tmp.hcat = i;

               /* calculate atom coordinates within large system */
               vec_t sys_coord;
               sys_coord.x = tmp.pos.x / double(mat.ucd.x) / double(sd.x);
               sys_coord.y = tmp.pos.y / double(mat.ucd.y) / double(sd.y);
               sys_coord.z = tmp.pos.z / double(mat.ucd.z) / double(sd.z);

               /* pin the centre slice in-plane */
               if (sys_coord.x > 0.49 && sys_coord.x < 0.51) {
                  tmp.element = "H";
                  n_pin ++;
                  tmp.mat += elements.size();
               }

               /* we want to split the system into two halves */
               if (sys_coord.x > 0.51) {
                  tmp.mat += 2*elements.size();
               }

               pinxyz << tmp.element << "\t"
                      << sys_coord.x << "\t"
                      << sys_coord.y << "\t"
                      << sys_coord.z << "\n";

               /* output to unit cell file */
               dwucf << tmp.gid << "\t"
                  << sys_coord.x << "\t"
                  << sys_coord.y << "\t"
                  << sys_coord.z << "\t"
                  << tmp.mat << "\t"
                  << 0 << "\t"
                  << tmp.hcat << "\n";

               dwsystem[i][j][k].push_back(tmp);
            }
         }
      }
   }

   std::cout << "number of pinned atoms: " << n_pin << ", " <<
                n_pin/float(natoms)*100 << "% of atoms" << std::endl;

   }


   dwucf << "# interactions n exctype, id i j dx dy dz Jij\n";

   /*******************************************************/
   /**  determine interactions for every atom in system  **/
   /**         using pre-calculated neighbour list       **/
   /*******************************************************/

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
                  if (dwsystem[i][j][k][atom].aid == uc_interactions[p].i.aid) {

                     pair_t tmp;

                     tmp.iid = int_counter;
                     tmp.i.gid = dwsystem[i][j][k][atom].gid;

                     tmp.i.mat = dwsystem[i][j][k][atom].mat;
                     tmp.i.element = dwsystem[i][j][k][atom].element;

                     tmp.j.element = uc_interactions[p].j.element;

                     tmp.exchange = uc_interactions[p].exchange;

                     /* assume atom j is within system to begin with */
                     tmp.ucd.x = 0;
                     tmp.ucd.y = 0;
                     tmp.ucd.z = 0;

                     /* check if atom j is within system boundaries */
                     int ucx = i + uc_interactions[p].ucd.x;
                     int ucy = j + uc_interactions[p].ucd.y;
                     int ucz = k + uc_interactions[p].ucd.z;

                     /* if any of these conditions are met
                      * then atom is out of bounds */

                     if (ucx < 0 || ucx >= sd.x ||
                         ucy < 0 || ucy >= sd.y ||
                         ucz < 0 || ucz >= sd.z ) {

                        /* periodic boundaries conditions */
                        if ( ucx < 0 ) {
                           ucx += sd.x;
                           tmp.ucd.x = -1;

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
                           tmp.ucd.y = -1;
                        }

                        if ( ucz < 0 ) {
                           tmp.exchange *= -1; /* antiferro periodic boundaries */
                           ucz += sd.z;
                           tmp.ucd.z = -1;
                        }

                        if ( ucx >= sd.x ) {
                           ucx -= sd.x;
                           tmp.ucd.x = 1;
                        }

                        if ( ucy >= sd.y ) {
                           ucy -= sd.y;
                           tmp.ucd.y = 1;
                        }

                        if ( ucz >= sd.z ) {
                           tmp.exchange *= -1; /* antiferro periodic boundaries */
                           ucz -= sd.z;
                           tmp.ucd.z = 1;
                        }

                        /* having changed uc coordinates obtain j.gid */
                        tmp.j.gid = dwsystem[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

                     }

                     /* else it is within bounds so simply extract j.gid */
                     else tmp.j.gid = dwsystem[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

                     uc_interactions.push_back(tmp);

                     /* increment interaction id */
                     int_counter ++;

                  }
               }
            }
         }
      }
   }

   // std::cout << "number of interactions in system: " << interactions.size() << "\n\n";

   dwucf << uc_interactions.size() << "\tisotropic\n";

   // output interaction info to file
   for (int i=0; i<uc_interactions.size(); i++)

      dwucf << uc_interactions[i].iid << "\t"
            << uc_interactions[i].i.gid << "\t"
            << uc_interactions[i].j.gid << "\t"
            << uc_interactions[i].ucd.x << "\t"
            << uc_interactions[i].ucd.y << "\t"
            << uc_interactions[i].ucd.z << "\t"
            << uc_interactions[i].exchange << "\n";

   std::cout << "number of interactions found: " << uc_interactions.size() << std::endl;
   std::cout << "interaction data output to file 'domainwall.ucf'\n";

   dwucf.close();

   return EXIT_SUCCESS;
}

int set_default_parameters() {

   /* cut off radii */
   sys.rcut_tt = 5.0;
   sys.rcut_rt = 5.0;

   /* exchange parameters */
   sys.tt_factor = 1.0;
   sys.rt_factor = 1.0;
   double rt_constant = 1.0;
   double tt_constant = 1.0;

   /* doping concentrations */
   sys.zr_concentration = 0.0;
   sys.ti = 0.0;

   /* flags */
   sys.tracking = false;
   sys.domainwall = false;
   sys.centrepin = false;

   /* domain wall parameters */
   sys.dw_dim.x = 0.0;
   sys.dw_dim.y = 0.0;
   sys.dw_dim.z = 0.0;

   return EXIT_SUCCESS;

}

// for now this function will include tracked cell calculation */
int generate_large_system(std::vector<pair_t>& uc_interactions,
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

         pair_t temp;

         temp.iid = int_counter;
         temp.i.gid = system[i][j][k][atom].gid;

         temp.i.mat = system[i][j][k][atom].mat;
         temp.i.element = system[i][j][k][atom].element;

         temp.j.element = uc_interactions[p].j.element;

         temp.exchange = uc_interactions[p].exchange;

         /* assume atom j is within system to begin with */
         temp.ucd.x = 0;
         temp.ucd.y = 0;
         temp.ucd.z = 0;

         /* check if atom j is within system boundaries */
         int ucx = i + uc_interactions[p].ucd.x;
         int ucy = j + uc_interactions[p].ucd.y;
         int ucz = k + uc_interactions[p].ucd.z;

         /* if any of these conditions satisfied
          * then atom is out of bounds
          */

         if (ucx < 0 || ucx >= sd.x ||
             ucy < 0 || ucy >= sd.y ||
             ucz < 0 || ucz >= sd.z  ) {

            /* periodic boundaries conditions */
            if ( ucx < 0 ) {
               ucx += sd.x;
               temp.ucd.x = -1;

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
                 temp.ucd.y = -1;
            }

            if ( ucz < 0 ) {
                 ucz += sd.z;
                 temp.ucd.z = -1;
            }

            if ( ucx >= sd.x ) {
                 ucx -= sd.x;
                 temp.ucd.x = 1;
            }

            if ( ucy >= sd.y ) {
                 ucy -= sd.y;
                 temp.ucd.y = 1;
            }

            if ( ucz >= sd.z ) {
                 ucz -= sd.z;
                 temp.ucd.z = 1;
            }

            /* having changed uc coordinates obtain j.gid */
            temp.j.gid = system[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

         }

         /* else it is within bounds so simply extract j.gid */
         else temp.j.gid = system[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

         uc_interactions.push_back(temp);

         /* increment interaction id */
         int_counter ++;

      }

      std::cout << "number of interactions in system: " << uc_interactions.size() << "\n\n";

      largeucf << uc_interactions.size() << "\tisotropic\n";

      // output interaction info to file
      for (int i=0; i<uc_interactions.size(); i++)

      largeucf << uc_interactions[i].iid << "\t"
          << uc_interactions[i].i.gid << "\t"
          << uc_interactions[i].j.gid << "\t"
          << uc_interactions[i].ucd.x << "\t"
          << uc_interactions[i].ucd.y << "\t"
          << uc_interactions[i].ucd.z << "\t"
          << uc_interactions[i].exchange << "\n";
   }

   largeucf.close();

   return EXIT_SUCCESS;
}
