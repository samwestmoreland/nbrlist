#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "../hdr/classes.hpp"

/* function prototypes */
double calculate_rij(vec_t& i, vec_t& j);     /* distance calculation */

double calculate_jij(std::string const& i_type,
                     std::string const& j_type,
                     double rij,
                     double tt_factor,
                     double rt_factor,
                     int exchange_fn);

double jij_ndfeb(std::string const& i_type,
                 std::string const& j_type,
                 double rij,
                 double tt_factor,
                 double rt_factor);

int initialise_material(int material, std::string const& material_string, double zr_content);
void array_to_rasmol(std::vector<atom_t> array, std::string const& arrayname);
int determine_material_id(std::string const& in_material);

int calculate_interactions(int exchange_fn,
                           double tt_factor,
                           double rt_factor,
                           double rcut);

vec_t get_uc_dimensions_from_zr_content(double zr_content);
int output_materials(std::vector<material_t>& materials);

int generate_large_system(std::vector<int_t>& uc_interactions,
                          std::vector<atom_t>& unitcell,
                          std::vector < std::vector < std::vector < std::vector <atom_t> > > >& system, vec_t ucd,
                          int n_tracked_cells,
                          int n_materials,
                          int system_dimension);

double ask_for_cutoff(double rcut);
double ask_for_rt(double rt_factor);
double ask_for_tt(double tt_factor);
int ask_for_system_dimensions(int system_dimension);
int ask_for_n_tracked_cells(int n_tracked_cells);
double ask_for_zr_content(double zr_content);
int ask_for_uc_config(int config);

std::string generate_filename(std::string const& material_string);

/* arrays (global) */
std::vector<atom_t> unitcell;
std::vector<atom_t> supercell;
std::vector<int_t> uc_interactions;
std::vector<material_t> materials;
vec_t ucd;

/* output file stream */
std::ofstream outfile;

int main (int argc, char *argv[]) {

   /* material options are:
      * ndfeb
      * ndfe12
      * bccfe
      * smfe12
      * smzrfe12
      * interface
      * interface_mirror
   */

   int material = atoi(argv[1]);

   /* global variables */
   double rcut = 0;
   bool tracking = false;

   double tt_factor = 0;
   double rt_factor = 0;

   int system_dimension = 0;
   int n_tracked_cells = 0;

   double zr_content = 0;
   int config = 0;

   int exchange_fn;

   std::string material_string;

   rcut = ask_for_cutoff(rcut);
   tt_factor = ask_for_tt(tt_factor);
   rt_factor = ask_for_rt(rt_factor);
   exchange_fn = material;

   if (tracking) {
      system_dimension = ask_for_system_dimensions(system_dimension);
      n_tracked_cells = ask_for_n_tracked_cells(n_tracked_cells);
   }

   switch(material) {

      case 1 :    // bccfe
         material_string = "bccfe";
         initialise_material(material, material_string, zr_content);
         break;

      case 2 :    // ndfeb
         material_string = "ndfeb";
         initialise_material(material, material_string, zr_content);
         break;

      case 3 :    // ndfe12
         material_string = "ndfe12";
         initialise_material(material, material_string, zr_content);
         break;

      case 4 :    // smfe12
         material_string = "smfe12";
         initialise_material(material, material_string, zr_content);
         break;

      case 5 :    // smzrfe12
         ask_for_zr_content(zr_content);
         ask_for_uc_config(config);
         material_string = "smzrfe12";
         initialise_material(material, material_string, zr_content);
         break;

      case 6 :    // interface
         material_string = "interface";
         initialise_material(material, material_string, zr_content);
         break;

      case 7 :    // interface_mirror
         material_string = "interface_mirror";
         initialise_material(material, material_string, zr_content);
         break;

      default :
         std::cout << "invalid option. exiting.\n";
         exit(EXIT_SUCCESS);
   }

   if (argc > 2)
   {
      std::string argv2 = argv[2];  // convert char* to string
      if (argv2 == "tracking")
         tracking = true;
   }

   /* *************************
    * **** read parameters ****
    * *************************/

   // if (material == "smzrfe12")
   // {
   //     std::cout << "input Zr parameters: zr_content (x), unit cell configuration\n";
   //     std::cin >> zr_content >> config;
   // }

   /* ************************
    * *** print parameters ***
    * ************************/

   std::cout << "\nuser inputted parameters\n";
   std::cout << "material: " << material_string << std::endl;
   std::cout << "cut-off radius: " << rcut << std::endl;

   // if (tracking)
   // {
   //     std::cout << "tracking system dimension: " << system_dimension << std::endl;
   //     std::cout << "no. of tracked cells: " << n_tracked_cells << std::endl;
   // }

   // if (material == "smzrfe12")
   // {
   //     std::cout << "Zr content (x): " << zr_content << std::endl;
   //     std::cout << "unit cell configuration: " << config << std::endl;
   // }

   std::cout << std::endl;

   array_to_rasmol(unitcell, "unitcell");

   std::cout << "number of materials: " << materials.size() << "\n\n\t";

   /* output names of materials */
   output_materials(materials);

   /* this function generates the supercell and fills the interactions array */
   calculate_interactions(exchange_fn, tt_factor, rt_factor, rcut);

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

          else std::cout << species_interactions[i][j] << "\t";
       }

       std::cout << std::endl;
   }

   std::cout << std::endl;

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

std::string generate_filename(std::string const& material_string) {

   std::string filename;

   filename = "./coordinates/" + material_string + ".coords";

   return filename;
}

double ask_for_zr_content(double zr_content) {

   std::cout << "enter zr content (x)\n";
   std::cin >> zr_content;

   return zr_content;
}

int ask_for_uc_config(int config) {

   std::cout << "enter unit cell configuration\n";
   std::cin >> config;

   return config;
}

int ask_for_system_dimensions(int system_dimension) {

   std::cout << "enter system dimension\n";
   std::cin >> system_dimension;

   return system_dimension;
}

int ask_for_n_tracked_cells(int n_tracked_cells) {

   std::cout << "enter number of tracked cells\n";
   std::cin >> n_tracked_cells;

   return n_tracked_cells;

}

double ask_for_rt(double rt_factor) {

   std::cout << "enter exchange: RE-TM\n";
   std::cin >> rt_factor;

   return rt_factor;
}

double ask_for_tt(double tt_factor) {

   std::cout << "enter exchange: TM-TM\n";
   std::cin >> tt_factor;

   return tt_factor;
}

double ask_for_cutoff(double rcut) {

   std::cout << "enter cut-off radius\n";
   std::cin >> rcut;

   return rcut;
}

/* for now this function will include tracked cell calculation */
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

   std::ofstream ucf ("large.ucf");

   /* output coordinates to unit cell file */
   ucf << "# Unit cell size:\n"
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

               atom_t tmp;
               vec_t uc;
               uc.x = i;
               uc.y = j;
               uc.z = k;

               tmp.aid = unitcell[atom].aid;
               tmp.gid = gid_counter;

               tmp.element = unitcell[atom].element;

               tmp.mat = unitcell[atom].mat;

               /*** cell tracking ***/


               if (  (i <  (sd.x+n_tracked_cells)*0.5)
                  && (i >= (sd.x-n_tracked_cells)*0.5)

                  && (j <  (sd.y+n_tracked_cells)*0.5)
                  && (j >= (sd.y-n_tracked_cells)*0.5)

                  && (k <  (sd.z+n_tracked_cells)*0.5)
                  && (k >= (sd.z-n_tracked_cells)*0.5) ) {

                  tmp.mat += n_materials;
                  tmp.element = 'h';
                  tracked ++;
               }

               tmp.pos = unitcell[atom].pos + uc*ucd;
               gid_counter ++;

               /* output coordinates to file */
               sysmol << tmp.element << "\t"
                      << tmp.pos.x << "\t"
                      << tmp.pos.y << "\t"
                      << tmp.pos.z << "\n";

               /* calculate atom coordinates within large system */
               vec_t sys_co;
               sys_co.x = tmp.pos.x / double(ucd.x) / double(sd.x);
               sys_co.y = tmp.pos.y / double(ucd.y) / double(sd.y);
               sys_co.z = tmp.pos.z / double(ucd.z) / double(sd.z);

               /* output to unit cell file */
               ucf << tmp.gid << "\t"
                   << sys_co.x << "\t"
                   << sys_co.y << "\t"
                   << sys_co.z << "\t"
                   << tmp.mat << "\t"
                   << 0 << "\t"
                   << 0 << "\n";

               system[i][j][k].push_back(tmp);

            }
         }
      }
   }

   std::cout << "total tracked atoms: " << tracked << "\n";
   std::cout << "total tracked cells (3D): " << tracked/unitcell.size() << "\n";

   ucf << "# interactions n exctype, id i j dx dy dz Jij\n";

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

         int_t tmp;

         tmp.iid = int_counter;
         tmp.i.gid = system[i][j][k][atom].gid;

         tmp.i.mat = system[i][j][k][atom].mat;
         tmp.i.element = system[i][j][k][atom].element;

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
          * then atom is out of bounds
          */

         if (ucx < 0 || ucx >= sd.x ||
             ucy < 0 || ucy >= sd.y ||
             ucz < 0 || ucz >= sd.z  ) {

            /* periodic boundaries conditions */
            if ( ucx < 0 ) {
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
            tmp.j.gid = system[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

         }

         /* else it is within bounds so simply extract j.gid */
         else tmp.j.gid = system[ucx][ucy][ucz][uc_interactions[p].j.aid].gid;

         interactions.push_back(tmp);

         /* increment interaction id */
         int_counter ++;

      }

      std::cout << "number of interactions in system: " << interactions.size() << "\n\n";

      ucf << interactions.size() << "\tisotropic\n";

      // output interaction info to file
      for (int i=0; i<interactions.size(); i++)

      ucf << interactions[i].iid << "\t"
          << interactions[i].i.gid << "\t"
          << interactions[i].j.gid << "\t"
          << interactions[i].disp.x << "\t"
          << interactions[i].disp.y << "\t"
          << interactions[i].disp.z << "\t"
          << interactions[i].exchange << "\n";
   }

   ucf.close();

   return EXIT_SUCCESS;
}

int output_materials(std::vector<material_t>& materials) {

    for (int i=0; i<materials.size(); ++i)
        std::cout << materials[i].name << "\t" << materials[i].id << "\n\t";

    return EXIT_SUCCESS;
}

int initialise_material(int material, std::string const& material_string, double zr_content) {

   std::cout << "initialising material \'" << material_string << "\'\n";

   switch(material) {

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
   std::string filename = generate_filename(material_string);

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

   std::cout << "reading in unit cell coordinates\n";

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
            temp_mat.id = determine_material_id(temp_mat.name);
            temp.mat = temp_mat.id;

            /* scale atom coordinates */
            temp.pos = temp.pos * ucd;

            unitcell.push_back(temp);
      }

   std::cout << "atoms read in: " << unitcell.size() << std::endl;

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

vec_t get_uc_dimensions_from_zr_content(double zr_content) {

    double atom_rad = (zr_content - 8.22624)/(-4.51952);
    double a = 0.0929088 * atom_rad + 0.830892;
    double c = -0.00593675 * zr_content + 1;

    a *= 8.497;
    c *= 4.687;

    vec_t ucd;
    ucd.x = a*2;
    ucd.y = a;
    ucd.z = c;
    return ucd;
}

/* returns a material */
int determine_material_id(std::string const& in_material) {

   /* define out_material which is a material_t to pushback into material array */
   material_t out_material;
   out_material.name = in_material;
   out_material.id = 0;

   /* if this is the first material */
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

int calculate_interactions(int exchange_fn, double tt_factor, double rt_factor, double rcut) {

    std::cout << "\npopulating super cell\n";

    populate_supercell();

    std::cout << "calculating interactions\n";

    std::cout << "using exchange function: " << exchange_fn << std::endl;

    if (exchange_fn == 2)
        std::cout <<
            "TM-TM exchange factor = " << tt_factor << "\n" <<
            "RE-TM exchange factor = " << rt_factor << "\n";

    /* central cell location */
    int start = (supercell.size()-unitcell.size())/2;
    int end   = (supercell.size()+unitcell.size())/2;

    int interaction_count = 0;

    /* loop through central cell */
    for (int i=start; i<end; ++i) {

        /* loop through super cell (looking for atom j) */
        for (int j=0; j<supercell.size(); ++j) {

            /* calculate interatomic distance */
            double rij = calculate_rij(supercell[i].pos, supercell[j].pos);

            /* if distance less than rcut and not same atom */
            if (rij < rcut && rij > 1e-30) {

                /* create an interaction */
                int_t temp;
                temp.i = supercell[i];
                temp.j = supercell[j];

                temp.iid = interaction_count;

                /* unitcell displacement */
                temp.disp = temp.j.uc - temp.i.uc;

                /* calculate exchange energy */
                temp.exchange = calculate_jij(temp.i.element, temp.j.element, rij, tt_factor, rt_factor, exchange_fn);

                /* put interaction into array */
                if (temp.exchange!=0) {

                    uc_interactions.push_back(temp);
                    interaction_count ++;

                }
            }
        }
    }

    std::cout
        << "interactions: "
        << uc_interactions.size() << std::endl;

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
        std::cout << filename << " file generated.\n";
    }
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
         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;
         double d = 2.0; /* factor 2 to correct for Hamiltonian */
         return d*c*(a/(rij*rij*rij)-b);
      }
         break;

      case 2 : // ndfeb

         return jij_ndfeb(i_type, j_type, rij, tt_factor, rt_factor);
         break;

      case 3 : { // ndfe12

         const double Fe_ratio_ndfe12=1.15*1.07692307692; // 560/520 = 1.07692307692
         const double J0Nd_ndfe12=Fe_ratio_ndfe12*4.06835e-20/16.0;

         /* nd-nd */
         if (i_type=="Nd" && j_type=="Nd") return 0.0;

         else return 0.0;
      }
         break;

      case 4 : { /* smfe12 */

         double a = 60.5033;
         double b = 0.862717;
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

      case 5 : /* smzrfe12 */

         return jij_ndfeb(i_type, j_type, rij, tt_factor, rt_factor);
         break;

      case 6 : /* interface */

         return jij_ndfeb(i_type, j_type, rij, tt_factor, rt_factor);
         break;


      case 7 : /* interface mirror */

         return jij_ndfeb(i_type, j_type, rij, tt_factor, rt_factor);
         break;

   } /* end of switch */

} /* end of calculate_jij function */

double jij_ndfeb(std::string const& i_type,
                 std::string const& j_type,
                 double rij,
                 double tt_factor,
                 double rt_factor) {

   /* richard's code */
   const double a = 36.9434;
   const double b = 1.25094;
   const double c = -0.229572;
   const double Fe_ratio_ndfeb=0.69618016759*1.07692307692; // 560/520 = 1.07692307692
   const double J0Nd_ndfeb=Fe_ratio_ndfeb*4.06835e-20/16.0;

   /* nd-nd */
   if (i_type == "Nd" && j_type == "Nd") return 0.0;

         // nd-fe (step function at r = 4A)
   else if ((i_type == "Fe" && j_type == "Nd") || (i_type == "Nd" && j_type == "Fe")) {
      if (rij <= 4.0) return 0.33 * J0Nd_ndfeb;       // ndfeb
      else return 0.0;
   }

   /* fe-fe (cutoff at r = 5.74A) */
   else if (i_type=="Fe" && j_type=="Fe") {
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

//    // ucd.x = 106.7;
//    // ucd.y = 36.595;
//    // ucd.z = 36.595;