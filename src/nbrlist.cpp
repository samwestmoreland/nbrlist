// scaling of m for generic bcc system

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "../hdr/classes.hpp"

/* function prototypes */
double calc_rij(vec_t i, vec_t j);     // distance calculation
double calc_jij(std::string i_type, std::string j_type, double rij);
int initialise_material(std::string material);
int initialise_data_structures();
void array_to_rasmol(std::vector<atom_t> array, std::string arrayname);
material_t determine_material_id(material_t material);
int calculate_interactions();

/* arrays */
std::vector<atom_t> unitcell;
std::vector<atom_t> supercell;
std::vector<int_t> uc_interactions;
std::vector<material_t> materials;
vec_t ucd;

/* output file stream */
std::ofstream outfile;

int main (int argc, char *argv[])
{
   /* material options are:
      * ndfeb
      * ndfe12
      * bccfe
      * smfe12
   */

   char * material = argv[1];
   std::string str(material);
   initialise_material(material);

   std::cout << "number of materials: " << materials.size() << std::endl;

   calculate_interactions();

return 0;

}

int initialise_material(std::string material)
{
   std::cout << "initialising material " << material << std::endl;
   std::string filename = "coordinates/" + material + ".coords";

   std::ifstream infile (filename.c_str());

   /* check if file opened */
   if (!(infile.is_open()))
   {
      std::cerr   << "couldn't open coordinate file" << std::endl;
      std::cerr   << "exiting" << std::endl;
      std::exit(EXIT_FAILURE);
   }

   if (material == "ndfeb")
   {
       ucd.x = 8.8;
       ucd.y = 8.8;
       ucd.z = 12.2;
   }

   else if (material == "bccfe")
   {
       ucd.x = 2.856;
       ucd.y = 2.856;
       ucd.z = 2.856;
   }

   else if (material == "test")
   {
       ucd.x = 1;
       ucd.y = 1;
       ucd.z = 1;
   }

   else
   {
       std::cout << "material not found\nexiting program\n";
       std::exit(EXIT_FAILURE);
   }

   std::cout
      << "unit cell dimensions: "
      << ucd.x << ", "
      << ucd.y << ", "
      << ucd.z << "\n";

   atom_t temp;
   int atom_count = 0;
   int material_count = 0;

   std::cout << "\nreading in unit cell coordinates\n";

   while (infile
      >> temp.element
      >> temp.pos.x
      >> temp.pos.y
      >> temp.pos.z)
      {
            /* assign atom id */
            temp.aid = atom_count;
            atom_count ++;

            /* determine material id */
            material_t temp_mat;
            temp_mat.name = temp.element;
            temp_mat = determine_material_id(temp_mat);
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

   for (int i=0; i<unitcell.size(); ++i)
   {
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

material_t determine_material_id(material_t in_material)
{
    material_t out_material;
    out_material.name = in_material.name;
    out_material.id = 0;

    if (materials.size() == 0)
    {
        materials.push_back(out_material);
        return out_material;
    }

    for (int i=0; i<materials.size(); ++i)
    {
        if (in_material.name == materials[i].name)
        {
            out_material.id = materials[i].id;
            return out_material;
        }
    }

    /* if we get to here, material has not been found */
    out_material.id = materials.size();
    materials.push_back(out_material);

    return out_material;
}

void populate_supercell()
{
    /* loop through dimensions */
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)

                /* loop through atoms in unitcell */
                for (int atom=0; atom<unitcell.size(); atom++)
                {
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

                    /* place atom in array */
                    supercell.push_back(temp);
                }

    std::cout
        << "atoms in super cell: "
        << supercell.size() << std::endl;

    //    array_to_rasmol(supercell, "supercell");

}

int calculate_interactions()
{
    std::cout << "\npopulating super cell\n";

    populate_supercell();

    /* cut-off radius in angstroms */
    double rcut = 5.0;

    /* central cell location */
    int start = (supercell.size()-unitcell.size())/2;
    int end   = (supercell.size()+unitcell.size())/2;

    int interaction_count = 0;

    /* loop through central cell */
    for (int i=start; i<end; ++i)
    {
        /* loop through super cell (looking for atom j) */
        for (int j=0; j<supercell.size(); ++j)
        {
            /* calculate interatomic distance */
            double rij = calc_rij(supercell[i].pos, supercell[j].pos);

            /* if distance less than rcut and not same atom */
            if (rij < rcut && rij > 1e-30)
            {
                /* create an interaction */
                int_t temp;
                temp.i = supercell[i];
                temp.j = supercell[j];

                temp.iid = interaction_count;

                /* unitcell displacement */
                temp.disp = temp.j.uc - temp.i.uc;

                /* calculate exchange energy */
                temp.exchange = calc_jij(temp.i.element, temp.j.element, rij);

                /* put interaction into array */
                if (temp.exchange!=0)
                {
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

void array_to_rasmol(std::vector<atom_t> array, std::string arrayname)
{
    std::string filename = arrayname + ".xyz";
    std::ofstream rasmol (filename.c_str());

    /* if file open */
    if (rasmol)
    {
        rasmol << array.size() << "\n\n";

        for (int i=0; i<array.size(); ++i)
            rasmol
                << array[i].element << "\t"
                << array[i].pos.x   << "\t"
                << array[i].pos.y   << "\t"
                << array[i].pos.z   << "\t"
                << array[i].uc.x    << "\t"
                << array[i].uc.y    << "\t"
                << array[i].uc.z    << "\n";

        rasmol.close();
        std::cout << filename << " file generated.\n";
    }
}

/* exchange function */
double calc_jij (std::string i_type, std::string j_type, double rij)
{

    double A=36.9434;
    double B=1.25094;
    double C=-0.229572;
    const double Fe_ratio_ndfeb=0.69618016759*1.07692307692; // 560/520 = 1.07692307692
    const double Fe_ratio_ndfe12=1.15*1.07692307692; // 560/520 = 1.07692307692
    // ^^ ndfeb

    const double J0Nd_ndfeb=Fe_ratio_ndfeb*4.06835e-20/16.0;
    const double J0Nd_ndfe12=Fe_ratio_ndfeb*4.06835e-20/16.0;

    // Nd-Nd
    if(i_type=="Nd" && j_type=="Nd"){
        return 0.0;
    }

    // Nd-Fe (step function at r = 4A)
    if((i_type=="Fe" && j_type=="Nd") || (i_type=="Nd" && j_type=="Fe"))
    {
        if(rij<=4.0) return 0.33*J0Nd_ndfeb;       // ndfeb
        else return 0.0;
    }

    // Fe-Fe (cutoff at r = 5.74A)
    else if(i_type=="Fe" && j_type=="Fe"){
        //if(rij<=5.0) return -2.0*2.179872e-21*(A*exp(-B*rij)+C);
        // Correct for Tc = 600
        if(rij<=5.0) return 2.0*2.179872e-21*(A*exp(-B*rij)+C)*Fe_ratio_ndfeb;
        else return 0.0;

        //        if (rij <= 2.6) return 3.117407e-21;
        //        else if (rij <= 3.0) return 1.7742e-21;
        //        else return 0.0;
    }

    // B-x
    else return 0.0;

}

// function to calculate distance
double calc_rij (vec_t i, vec_t j)
{
    vec_t d = j-i;
    double distance = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);

    return distance;
}

//    // ucd.x = 106.7;
//    // ucd.y = 36.595;
//    // ucd.z = 36.595;

//    /*
//    * generate large system
//    */
//
//    // system dimensions
//    vec_t sd;
//    sd.x = 1;
//    sd.y = 1;
//    sd.z = 1;
//
//    std::cout << "\ndimensions of large system [A]: ("
//    << sd.x * ucd.x << ", "
//    << sd.y * ucd.y << ", "
//    << sd.z * ucd.z << ")\n";
//
//    int ns = sd.x * sd.y * sd.z * unitcell.size(); // number of atoms in system
//    std::cout << "number of atoms in system: " << ns << "\n";
//
//    // open file for rasmol output of system
//    std::ofstream sysxyz ("system.xyz");
//    sysxyz << ns << "\n\n";
//
//    // open file for ucf output
//    std::ofstream ucf ("ndfeb.ucf");            // change
//
//    ucf << "# Unit cell size:\n"
//    <<  sd.x*ucd.x << "\t"
//    <<  sd.y*ucd.y << "\t"
//    <<  sd.z*ucd.z << "\n"
//    << "# Unit cell vectors:\n"
//    << "1.0  0.0  0.0\n"
//    << "0.0  1.0  0.0\n"
//    << "0.0  0.0  1.0\n"
//    << "# Atoms num, id cx cy cz mat lc hc\n"
//    << ns << "\n";
//
//    // array to hold all atoms in system
//    std::vector < std::vector < std::vector < std::vector <atom_t> > > > system;
//
//    // global id counter
//    int gid_counter = 0;
//    int tracked = 0;
//
//    /* resize vectors */
//    system.resize(sd.x);
//    for (int i=0; i<sd.x; i++) {
//       system[i].resize(sd.y);
//       for (int j=0; j<sd.y; j++) {
//          system[i][j].resize(sd.z);
//          for (int k=0; k<sd.z; k++) {
//
//             /* loop through unitcell atoms */
//             for (int atom=0; atom<unitcell.size(); ++atom) {
//
//                atom_t tmp;
//                vec_t uc;
//                uc.x = i;
//                uc.y = j;
//                uc.z = k;
//
//                tmp.aid = unitcell[atom].aid;
//                tmp.gid = gid_counter;
//
//                tmp.element = unitcell[atom].element;
//
//                tmp.mat = unitcell[atom].mat;
//
//                /*
//                * cell tracking
//                */
//
//                // if ( (i<((sd.x+ntr/2)-sd.x/2))
//                //    &&(i>((sd.x-ntr/2)-sd.x/2))
//                //    &&(j<((sd.y+ntr/2)-sd.y/2))
//                //    &&(j>((sd.y-ntr/2)-sd.y/2))
//                //    &&(k<((sd.z+ntr/2)-sd.z/2))
//                //    &&(k>((sd.z-ntr/2)-sd.z/2)))
//                // {
//                //     tmp.mat = 1;
//                //     tmp.element = 'h';
//                //     ++tracked;
//                // }
//
//                tmp.pos = unitcell[atom].pos + uc*ucd;
//                ++gid_counter;
//
//                sysxyz << tmp.element << "\t"
//                << tmp.pos.x << "\t"
//                << tmp.pos.y << "\t"
//                << tmp.pos.z << "\n";
//
//                // normalise coordinates to large system dimensions
//                double sysx = tmp.pos.x / double(ucd.x) / double(sd.x);
//                double sysy = tmp.pos.y / double(ucd.y) / double(sd.y);
//                double sysz = tmp.pos.z / double(ucd.z) / double(sd.z);
//
//                // height categorisation
//                double hinc = 6;            // change this
//                tmp.hcat = sysz*hinc;
//
//                // output to file
//                ucf << tmp.gid << "\t"
//                << sysx << "\t"
//                << sysy << "\t"
//                << sysz << "\t"
//                << tmp.mat << "\t"
//                << 0 << "\t"
//                << tmp.hcat << "\n";
//
//                system[i][j][k].push_back(tmp);
//
//             }
//          }
//       }
//    }
//
//    std::cout << "number of tracked cells: " << tracked/unitcell.size() << "\n";
//
//    sysxyz.close();
//
//    ucf << "# Interactions n exctype, id i j dx dy dz Jij\n";
//
//    // ===============================================
//    // determine interactions for every atom in system
//    // using pre-calculated neighbour list
//    // ===============================================
//
//    // initialise array to hold interactions for whole system
//    std::vector<int_t> interactions;
//
//    // interaction counter
//    int int_counter = 0;
//
//    // loop through every atom in system
//    for (int i=0; i<sd.x; i++)
//    for (int j=0; j<sd.y; j++)
//    for (int k=0; k<sd.z; k++)
//
//    // first atom
//    for (int atom=0; atom<unitcell.size(); atom++)
//
//    // loop through interaction information
//    for (int p=0; p<ucints.size(); p++) {
//
//       // if interaction info refers to correct atom
//       if (system[i][j][k][atom].aid == ucints[p].i.aid)
//       {
//          int_t tmp;
//
//          tmp.iid = int_counter;
//          tmp.i.gid = system[i][j][k][atom].gid;
//
//          tmp.i.mat = system[i][j][k][atom].mat;
//          tmp.i.element = system[i][j][k][atom].element;
//
//          tmp.j.element = ucints[p].j.element;
//
//          tmp.exchange = ucints[p].exchange;
//
//          // assume atom j is within system
//          // to begin with
//          tmp.disp.x = 0;
//          tmp.disp.y = 0;
//          tmp.disp.z = 0;
//
//          // check if atom j is within system boundaries
//          int ucx = i + ucints[p].disp.x;
//          int ucy = j + ucints[p].disp.y;
//          int ucz = k + ucints[p].disp.z;
//
//          // only one of these conditions
//          // needs to be satisfied for the
//          // atom to be out of bounds
//          if ( ucx < 0 || ucx >= sd.x ||
//             ucy < 0 || ucy >= sd.y ||
//             ucz < 0 || ucz >= sd.z )
//             {
//
//                // periodic boundaries conditions
//                if (ucx<0) {
//                   ucx += sd.x;
//                   tmp.disp.x = -1;
//
//                   //    // determine mat of atom j
//                   //    if (tmp.i.mat==5)
//                   //    {
//                   //        if (tmp.j.element=="Nd") tmp.j.mat = 5;
//                   //        else tmp.j.mat = 6;
//                   //    }
//                   //
//                   //    else if (tmp.i.mat==6)
//                   //    {
//                   //        if (tmp.j.element=="Fe") tmp.j.mat = 6;
//                   //        else tmp.j.mat = 5;
//                   //    }
//
//                }
//
//                if (ucy<0) {
//                   ucy += sd.y;
//                   tmp.disp.y = -1;
//                }
//
//                if (ucz<0) {
//                   ucz += sd.z;
//                   tmp.disp.z = -1;
//                }
//                if (ucx>=sd.x) {
//                   ucx -= sd.x;
//                   tmp.disp.x = 1;
//                }
//                if (ucy>=sd.y) {
//                   ucy -= sd.y;
//                   tmp.disp.y = 1;
//                }
//                if (ucz>=sd.z) {
//                   ucz -= sd.z;
//                   tmp.disp.z = 1;
//                }
//
//                // having changed uc coordinates
//                // obtain j.gid
//                tmp.j.gid = system[ucx][ucy][ucz][ucints[p].j.aid].gid;
//
//             }
//
//             // else it is within bounds
//             // so simply extract j.gid
//             else tmp.j.gid = system[ucx][ucy][ucz][ucints[p].j.aid].gid;
//
//             /*
//             // if at top afm boundary
//             // switch sign of exchange
//             // and amplify exchange to pin wall
//             if (k==sd[2])
//             {
//             if (tmp.dz == -1) tmp.exchange *= -100;
//          }
//          else if (k==sd[2]-1)
//          {
//          if (tmp.dz == 1) tmp.exchange *= -100;
//       }
//       */
//
//       interactions.push_back(tmp);
//
//       // increment interaction id
//       ++int_counter;
//
//    }
// }
//
// // std::cout << "number of interactions in system: " << interactions.size() << "\n\n";
//
// ucf << interactions.size() << "\t0\n";
//
// // output interaction info to file
// for (int i=0; i<interactions.size(); i++)
//
// ucf << interactions[i].iid << "\t"
// << interactions[i].i.gid << "\t"
// << interactions[i].j.gid << "\t"
// << interactions[i].disp.x << "\t"
// << interactions[i].disp.y << "\t"
// << interactions[i].disp.z << "\t"
// << interactions[i].exchange << "\n";
//
// ucf.close();
//
